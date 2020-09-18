/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

template<class CloudType>
Foam::label Foam::ParticleStatistic<CloudType>::maxStoredParcels() const
{
    return maxStoredParcels_;
}

template<class CloudType>
inline const Foam::Switch&
Foam::ParticleStatistic<CloudType>::resetOnWrite() const
{
    return resetOnWrite_;
}

template<class CloudType>
Foam::vector Foam::ParticleStatistic<CloudType>::areaTri(const point& a, const point& b, const point& c) const
{
    // Reference: 
    // OpenFOAM-6/src/OpenFOAM/meshes/primitiveShapes/triangle/triangleI.H function area()
    // 2019/10/25 20:15:14  zt
    return 0.5*((b-a)^(c-a));
}

template<class CloudType>
Foam::vector Foam::ParticleStatistic<CloudType>::areaFace(const Field<point>& p) const
{   
    // Reference:
    // OpenFOAM-6/src/OpenFOAM/meshes/meshShapes/face/face.C function area(const pointField& p)
    // 2019/10/25 20:15:18  zt
    const label nPoints = p.size();

    // Calculate the area by summing the face triangle areas.
    // Changed to deal with small concavity by using a central decomposition
    //

    // If the face is a triangle, do a direct calculation to avoid round-off
    // error-related problems
    //
    if (nPoints == 3)
    {
        return areaTri
            (
             p[0],
             p[1],
             p[2]
            );
    }

    label pI;

    point centrePoint = vector::zero;
    for (pI = 0; pI < nPoints; ++pI)
    {
        centrePoint += p[pI];
    }
    centrePoint /= nPoints;

    vector a = vector::zero;

    point nextPoint = centrePoint;

    for (pI = 0; pI < nPoints; ++pI)
    {
        if (pI < nPoints - 1)
        {
            nextPoint = p[pI + 1];
        }
        else
        {
            nextPoint = p[0];
        }

        // Note: for best accuracy, centre point always comes last
        //
        a += areaTri(p[pI], nextPoint, centrePoint);
    }

    return a;
}

// ************************************************************************* //
