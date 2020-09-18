/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ParticleStatistic.H"
#include "Pstream.H"
#include "stringListOps.H"
#include "ListOps.H"
#include "ListListOps.H"

// Added from ParticleCollector
#include "surfaceWriter.H"
#include "unitConversion.H"
#include "Random.H"
#include "triangle.H"
#include "cloud.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleStatistic<CloudType>::initPolygons
(
    const List<Field<point> >& polygons
)
{
    mode_ = mtPolygon;

    label nPoints = 0;
    forAll(polygons, polyI)
    {
        label np = polygons[polyI].size();
        if (np < 3)
        {
            FatalIOErrorIn
            (
                "Foam::ParticleStatistic<CloudType>::initPolygons()",
                this->coeffDict()
            )
                << "polygons must consist of at least 3 points"
                << exit(FatalIOError);
        }

        nPoints += np;
    }

    label pointOffset = 0;
    points_.setSize(nPoints);
    faces_.setSize(polygons.size());
    faceTris_.setSize(polygons.size());
    forAll(faces_, faceI)
    {
        const Field<point>& polyPoints = polygons[faceI];
        face f(identity(polyPoints.size()) + pointOffset);
        UIndirectList<point>(points_, f) = polyPoints;

        DynamicList<face> tris;
        f.triangles(points_, tris);
        faceTris_[faceI].transfer(tris);

        faces_[faceI].transfer(f);

        pointOffset += polyPoints.size();
    }
}
template<class CloudType>
void Foam::ParticleStatistic<CloudType>::initConcentricCircles()
{
    mode_ = mtConcentricCircle;

    vector origin(this->coeffDict().lookup("origin"));

    radius_ = this->coeffDict().lookup("radius");
    nSector_ = readLabel(this->coeffDict().lookup("nSector"));

    label nS = nSector_;

    vector refDir;
    if (nSector_ > 1)
    {
        refDir = this->coeffDict().lookup("refDir");
        refDir -= normal_[0]*(normal_[0] & refDir);
        refDir /= mag(refDir);
    }
    else
    {
        // set 4 quadrants for single sector cases
        nS = 4;

        vector tangent = vector::zero;
        scalar magTangent = 0.0;

        Random rnd(1234);
        while (magTangent < SMALL)
        {
            vector v = rnd.vector01();

            tangent = v - (v & normal_[0])*normal_[0];
            magTangent = mag(tangent);
        }

        refDir = tangent/magTangent;
    }

    scalar dTheta = 5.0;
    scalar dThetaSector = 360.0/scalar(nS);
    label intervalPerSector = max(1, ceil(dThetaSector/dTheta));
    dTheta = dThetaSector/scalar(intervalPerSector);

    label nPointPerSector = intervalPerSector + 1;

    label nPointPerRadius = nS*(nPointPerSector - 1);
    label nPoint = radius_.size()*nPointPerRadius;
    label nFace = radius_.size()*nS;

    // add origin
    nPoint++;

    points_.setSize(nPoint);
    faces_.setSize(nFace);

    coordSys_ = cylindricalCS("coordSys", origin, normal_[0], refDir, false);

    List<label> ptIDs(identity(nPointPerRadius));

    points_[0] = origin;

    // points
    forAll(radius_, radI)
    {
        label pointOffset = radI*nPointPerRadius + 1;

        for (label i = 0; i < nPointPerRadius; i++)
        {
            label pI = i + pointOffset;
            point pCyl(radius_[radI], degToRad(i*dTheta), 0.0);
            points_[pI] = coordSys_.globalPosition(pCyl);
        }
    }

    // faces
    DynamicList<label> facePts(2*nPointPerSector);
    forAll(radius_, radI)
    {
        if (radI == 0)
        {
            for (label secI = 0; secI < nS; secI++)
            {
                facePts.clear();

                // append origin point
                facePts.append(0);

                for (label ptI = 0; ptI < nPointPerSector; ptI++)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = ptIDs.fcIndex(i - 1) + 1;
                    facePts.append(id);
                }

                label faceI = secI + radI*nS;

                faces_[faceI] = face(facePts);
            }
        }
        else
        {
            for (label secI = 0; secI < nS; secI++)
            {
                facePts.clear();

                label offset = (radI - 1)*nPointPerRadius + 1;

                for (label ptI = 0; ptI < nPointPerSector; ptI++)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = offset + ptIDs.fcIndex(i - 1);
                    facePts.append(id);
                }
                for (label ptI = nPointPerSector-1; ptI >= 0; ptI--)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = offset + nPointPerRadius + ptIDs.fcIndex(i - 1);
                    facePts.append(id);
                }

                label faceI = secI + radI*nS;

                faces_[faceI] = face(facePts);
            }
        }
    }
}

template<class CloudType>
void Foam::ParticleStatistic<CloudType>::collectParcelPolygon
(
    const point& p1,
    const point& p2
) const
{
    forAll(faces_, facei)
    {
        const label facePoint0 = faces_[facei][0];

        const point& pf = points_[facePoint0];

        const scalar d1 = normal_[facei] & (p1 - pf);
        const scalar d2 = normal_[facei] & (p2 - pf);

        if (sign(d1) == sign(d2))
        {
            // Did not cross polygon plane
            continue;
        }

        // Intersection point
        const point pIntersect = p1 + (d1/(d1 - d2))*(p2 - p1);

        // Identify if point is within the bounds of the face. Create triangles
        // between the intersection point and each edge of the face. If all the
        // triangle normals point in the same direction as the face normal, then
        // the particle is within the face. Note that testing for pointHits on
        // the face's decomposed triangles does not work due to ambiguity along
        // the diagonals.
        const face& f = faces_[facei];
        const vector a = areaFace(points_);
        bool inside = true;
        for (label i = 0; i < f.size(); ++ i)
        {
            const label j = f.fcIndex(i);
            const triPointRef t(pIntersect, points_[f[i]], points_[f[j]]);
            if ((a & areaTri(pIntersect, points_[f[i]], points_[f[j]])) < 0)
            {
                inside = false;
                break;
            }
        }

        // Add to the list of hits
        if (inside)
        {
            hitFaceIDs_.append(facei);
        }
    }
}

template<class CloudType>
void Foam::ParticleStatistic<CloudType>::collectParcelConcentricCircles
(
    const point& p1,
    const point& p2
) const
{
    label secI = -1;

    const scalar d1 = normal_[0] & (p1 - coordSys_.origin());
    const scalar d2 = normal_[0] & (p2 - coordSys_.origin());
    
    if (sign(d1) == sign(d2))
    {
        // did not cross plane
        return;
    }

    // intersection point in cylindrical co-ordinate system
    const point pCyl = coordSys_.localPosition(p1 + (d1/(d1 - d2))*(p2 - p1));

    scalar r = pCyl[0];

    if (r < radius_.last())
    {
        label radI = 0;
        while (r > radius_[radI])
        {
            radI++;
        }

        if (nSector_ == 1)
        {
            secI = 4*radI;
        }
        else
        {
            scalar theta = pCyl[1] + constant::mathematical::pi;

            secI =
                nSector_*radI
              + floor
                (
                    scalar(nSector_)*theta/constant::mathematical::twoPi
                );
        }
    }

    hitFaceIDs_.append(secI);
}
// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleStatistic<CloudType>::write()
{
    List<List<string> > procData(Pstream::nProcs());
    procData[Pstream::myProcNo()] = patchData_;
    Pstream::gatherList(procData);

    if (Pstream::master())
    {
        const fvMesh& mesh = this->owner().mesh();

        // Create directory if it doesn't exist
        mkDir(this->outputTimeDir());

        OFstream patchOutFile
            (
             this->outputTimeDir()/"statistic" + ".dat",
             IOstream::ASCII,
             IOstream::currentVersion,
             mesh.time().writeCompression()
            );

        List<string> globalData;
        globalData = ListListOps::combine<List<string> >
            (
             procData,
             accessOp<List<string> >()
            );

        string header("# "+parcelType::propertyList_);
        patchOutFile<< header.c_str() << nl;

        forAll(globalData, i)
        {
            patchOutFile
                << globalData[i].c_str()
                << nl;
        }
    }

    patchData_.clearStorage();

}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleStatistic<CloudType>::ParticleStatistic
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    mode_(mtUnknown),
    parcelType_(this->coeffDict().lookupOrDefault("parcelType", -1)),
    points_(),
    faces_(),
    faceTris_(),
    nSector_(0),
    radius_(),
    coordSys_(false),
    normal_(),
    surfaceFormat_(this->coeffDict().lookup("surfaceFormat")),
    resetOnWrite_(this->coeffDict().lookup("resetOnWrite")),
    outputFilePtr_(),
    hitFaceIDs_(), 
    maxStoredParcels_(readScalar(this->coeffDict().lookup("maxStoredParcels")))
{
    normal_ /= mag(normal_);

    word mode(this->coeffDict().lookup("mode"));
    if (mode == "polygon")
    {
        List<Field<point> > polygons(this->coeffDict().lookup("polygons"));

        initPolygons(polygons);

        vector n0(this->coeffDict().lookup("normal"));
        normal_ = vectorField(faces_.size(), n0);
    }
    else if (mode == "polygonWithNormal")
    {
        List<Tuple2<Field<point>, vector> > polygonAndNormal
        (
            this->coeffDict().lookup("polygons")
        );

        List<Field<point> > polygons(polygonAndNormal.size());
        normal_.setSize(polygonAndNormal.size());

        forAll(polygons, polyI)
        {
            polygons[polyI] = polygonAndNormal[polyI].first();
            normal_[polyI] = polygonAndNormal[polyI].second();
            normal_[polyI] /= mag(normal_[polyI]) + ROOTVSMALL;
        }

        initPolygons(polygons);
    }
    else if (mode == "concentricCircle")
    {
        vector n0(this->coeffDict().lookup("normal"));
        normal_ = vectorField(1, n0);

        initConcentricCircles();
    }
    else
    {
        FatalIOErrorIn
        (
            "Foam::ParticleStatistic<CloudType>::ParticleStatistic"
            "("
                "const dictionary&,"
                "CloudType&, "
                "const word&"
            ")",
            this->coeffDict()
        )
            << "Unknown mode " << mode << ".  Available options are "
            << "polygon, polygonWithNormal and concentricCircle"
            << exit(FatalIOError);
    }

}


template<class CloudType>
Foam::ParticleStatistic<CloudType>::ParticleStatistic
(
    const ParticleStatistic<CloudType>& ps
)
:
    CloudFunctionObject<CloudType>(ps),
    mode_(ps.mode_),
    parcelType_(ps.parcelType_),
    points_(ps.points_),
    faces_(ps.faces_),
    faceTris_(ps.faceTris_),
    nSector_(ps.nSector_),
    radius_(ps.radius_),
    coordSys_(ps.coordSys_),
    normal_(ps.normal_),
    surfaceFormat_(ps.surfaceFormat_),
    resetOnWrite_(ps.resetOnWrite_),
    outputFilePtr_(),
    hitFaceIDs_(),
    maxStoredParcels_(ps.maxStoredParcels_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleStatistic<CloudType>::~ParticleStatistic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleStatistic<CloudType>::postMove
(
    parcelType& p,
    const label cellI,
    const scalar dt,
    const point& position0,
    bool& keepParticle
)
{
    if ((parcelType_ != -1) && (parcelType_ != p.typeId()))
    {
        return;
    }

    hitFaceIDs_.clear();

    switch (mode_)
    {
        case mtPolygon:
        {
            collectParcelPolygon(position0, p.position());
            break;
        }
        case mtConcentricCircle:
        {
            collectParcelConcentricCircles(position0, p.position());
            break;
        }
        default:
        {
        }
    }
    forAll(hitFaceIDs_, i)
    {
        //Info<<"We are in hitFaceId: "<<hitFaceIDs_.size()<<"\t"<<hitFaceIDs_<<"\t"<<mode_<<"\t"<<p.origId()
        //<<position0<<'\t'<<p.position()<<endl;
        if (patchData_.size() < maxStoredParcels_)
        {
            OStringStream data;
            data<< p;
            patchData_.append(data.str());
        }
    }
}
// ************************************************************************* //
