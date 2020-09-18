/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::ParticleStatistic

Description
    Standard run-time Particle Statistic
    For particles pass through the sampling plane, capture it.
    Saved data includes all particle properties:
      - position (x,y,z)
      - nParticle (n particles in each parcel)
      - d (diameter)
      - U (velocity)
      - rho (density)
      - T (Temperature)

SourceFiles
    ParticleStatistic.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleStatistic_H
#define ParticleStatistic_H

#include "CloudFunctionObject.H"
#include "cylindricalCS.H"
#include "face.H"
#include "Switch.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ParticleStatistic Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleStatistic
:
    public CloudFunctionObject<CloudType>
{
public:

    enum modeType
    {
        mtPolygon,
        mtConcentricCircle,
        mtUnknown
    };


private:

    // Private data

        typedef typename CloudType::particleType parcelType;

        //- Collector mode type
        modeType mode_;

        //- Index of parcel types to collect (-1 by default = all particles)
        const label parcelType_;

        //- List of points
        Field<point> points_;

        //- List of faces
        List<face> faces_;

        // Polygon collector

            //- Triangulation of faces
            List<List<face> > faceTris_;

        // Concentric circles collector

            //- Number of sectors per circle
            label nSector_;

            //- List of radii
            List<scalar> radius_;

            //- Cylindrical co-ordinate system
            cylindricalCS coordSys_;
        
        //- Polygon normal vector per face
        Field<vector> normal_;

        //- Surface output format
        const word surfaceFormat_;

        //- Flag to indicate whether data should be reset/cleared on writing
        Switch resetOnWrite_;

        //- Output file pointer
        autoPtr<OFstream> outputFilePtr_;

        //- Work list to store which faces are hit
        mutable DynamicList<label> hitFaceIDs_;

        //- Maximum number of parcels to store - set as a scalar for I/O
        scalar maxStoredParcels_;

        //- List of output data per patch
        /*List<DynamicList<string> > patchData_;*/
        DynamicList<string> patchData_;

    // Private Member Functions
        //- Initialise polygon collectors
        void initPolygons(const List<Field<point> >& polygons);

        //- Initialise concentric circle collectors
        void initConcentricCircles();

        //- Collect parcels in polygon collectors
        void collectParcelPolygon
        (
            const point& p1,
            const point& p2
        ) const;

        //- Collect parcels in concentric circle collectors
        void collectParcelConcentricCircles
        (
            const point& p1,
            const point& p2
        ) const;

protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("particleStatistic");


    // Constructors

        //- Construct from dictionary
        ParticleStatistic
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleStatistic(const ParticleStatistic<CloudType>& ps);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType> > clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType> >
            (
                new ParticleStatistic<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleStatistic();


    // Member Functions

        // Access

            //- Return maximum number of parcels to store per patch
            inline label maxStoredParcels() const;

            //- Return const access to the reset on write flag
            inline const Switch& resetOnWrite() const;

            inline vector areaTri(const point& a, const point& b, const point& c) const;
            
            inline vector areaFace(const Field<point>& p) const;

        // Evaluation
        
            //- Post-move hook
            virtual void postMove
            (
                parcelType& p,
                const label cellI,
                const scalar dt,
                const point& position0,
                bool& keepParticle
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleStatisticI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ParticleStatistic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
