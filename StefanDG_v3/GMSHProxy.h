#pragma once
#include <cstdint>
#include "Coordinates.h"
#include "LocalCoordinates3D.h"
#include "GeometryConstants.h"

namespace GMSHProxy
{
    const size_t TASK_DEFAULT = 0;
    const size_t N_TASK_DEFAULT = 1;

    void initialize(int argc = 0, char** argv = 0,
                    const bool readConfigFiles = true,
                    const bool run = false);
    
    void finalize();

    void open(const char modelName[]);

    void fltkRun();

    void free(void* p);

    const int FALSE_C = 0;
    const int TRUE_C = 1;

    const double DEFAULT_MESH_SIZE = 2.0;

    namespace model
    {
        const int DEFAULT_TAG = -1;
        const int UNSUPPORTED_TAG = 0;
        const int DEFAULT_DIM = -1;

        void getEntitiesForPhysicalName(const char name[], int* &dimTags, size_t &nDimTags);
        void getEntitiesForPhysicalName(const char name[], int*& dimTags, unsigned int& nEntities);
        void getPhysicalGroups(int* &dimTags, size_t &nDimTags, const int dim = DEFAULT_DIM);
        void getPhysicalGroups(int*& dimTags, unsigned int& nGroups, int dim = DEFAULT_DIM);
        void getPhysicalGroups(const int dim, int** tags, size_t* nTags);
        void getPhysicalName(const int dim, const int tag, char* &physicalName);
        void getEntitiesForPhysicalGroup(const int dim, const int physicalTag, int* &tags, size_t &nEntities);
        void getEntitiesForPhysicalGroup(const int dim, const int physicalTag, int*& entitiesTags, unsigned int& nEntities);
        void getEntities(int* &dimTags, size_t &nEntities, const int dim = DEFAULT_DIM);
        void getEntities(int*& dimTags, unsigned int& nEntities, const int dim = DEFAULT_DIM);
        void getEntities(const int dim, int** tags, size_t* nEntities);
        void getSurfaceUpEntities(const int surfaceTag, int*& upEntitiesTags, size_t& nUpEntities);
        void getPhysicalGroupsForEntity(const int dim, const int entityTag, int* &physicalTags, size_t &nPhysicalGroups);
        void getBoundaries(const int entityDim, const int entityTag, int*& boundariesDimTags, size_t& nBoundariesDimTags);
        void getBoundaries(const int entityDim, const int entityTag, int*& boundariesDimTags, unsigned int& nBoundaries);
        void getBoundariesFor3DEntity(const int entityTag, int*& boundariesTags, size_t& nBoundariesTags);
        void getSurfaceBoundariesNodes(const int entityTag, int*& nodesTags, size_t& nNodes);
        void getParametrizationOnSurface(const int surfaceTag, const double* coordinates, const size_t nCoordinates, double*& parametricCoordinates, size_t& nParametricCoordinates);
        void getNormalsByParametric(const int surfaceTag, const double* parametricCoordinates, const size_t nParametricCoords, double*& normalsCoordinates, size_t& nNormalsCoordinates);
        void getNormalsByCoordinates(const int surfaceTag, const double* coordinates, const size_t nCoordinates, double*& normalsCoordinates);

        int addDescreteEntity(int dim);
        int addDescreteEntity2D();
        void addNodes(const double* nodesCoordinateIt, const size_t nNodes, int* nodesTagsIt);

        void removeEntities(const int dim, const int* tagsIt, const size_t nTags);
        void getSurafceType(const int surfaceTag, char*& surfaceType);
        bool isSurfacePlane(const int surfaceTag);
        void get2DPhysicalNames(const int* physicalTagsIt, const size_t nTags, char** physicalNamesIt);
        void deletePhysicalNames(char** physicalNamesIt, const size_t nPhysicalNames);
        void getCurrent(char*& modelName);
        void setCurrent(const char modleName[]);
        void remove();
        void removeEntities(const int* dimTags, const size_t nDimTags, bool recursive = false);
        void remove2DEntities(const int* tags, const size_t nTags, bool recursive = false);

        void extractTags(const int* dimTagIt, const size_t nTags, int* tagIt);

        namespace mesh
        {
            const int TRIANGLE_TYPE = 2;
            const int TETRAHEDRON_TYPE = 4;

            const int TRIANGLE_FACE_TYPE = 3;

            void generate(int dim);
            void createFaces();
            void createFaces(const int* dimTags, const size_t nDimTags);

            void getJacobians(const int elementType, double** jacobians, double** determinants, double** coord, size_t* nElements, const int tag = DEFAULT_TAG);
            void getJacobians(const int elementType, double*& jacobians, double*& determinants, double*& coord, size_t& nElements, const int tag = DEFAULT_TAG);
            void getJacobians(const int elementType,
                double (*&transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                double*& determinants,
                Coordinates*& initPoints,
                size_t& nElements,
                const int tag = DEFAULT_TAG);

            void getInitJacobian(const int elementTag, double*& jacobians, double& determinant);

            void getTriangleJacobians(double (*&transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                double*& determinants,
                Coordinates*& initPoints,
                const int tag = DEFAULT_TAG);

            void getTetrahedronsJacobians(double (*&transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                double*& determinants,
                Coordinates*& initPoints,
                size_t& nElements,
                const int tag = DEFAULT_TAG);

            void getTetrahedronsJacobians(double(*&transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                                         double*& determinants,
                                         Coordinates*& initPoints,
                                         const int tag);

            void getTrianglesNodes(size_t(*&traianglesNodesTags)[constants::triangle::N_NODES], Coordinates(*&traianglesNodes)[constants::triangle::N_NODES], const int tag = DEFAULT_TAG);


            void getIntegrationPoints(const int elementType,
                const char integrationType[],
                double** localCoord,
                double** weights,
                uint8_t* nIntegrationSteps);

            void getElementsByType(const int elementType,
                size_t* &elementTags,
                size_t &nElement,
                size_t* &nodeTags,
                const int tag = DEFAULT_TAG,
                const size_t task = TASK_DEFAULT,
                const size_t nTasks = N_TASK_DEFAULT);

            void getTetrahedrons(size_t*& tetrahedronsTags,
                size_t& nTetrahedrons,
                size_t(*&tetrahedronsNodesTags)[constants::tetrahedron::N_NODES],
                const int tag = DEFAULT_TAG,
                const size_t task = TASK_DEFAULT,
                const size_t nTasks = N_TASK_DEFAULT);


            void getFacesByElements(const int elementType, const int faceType, size_t** nodesTags, size_t* nNodesTags, const int entityTag = DEFAULT_TAG);
            void getTetrahedronsFacesByElements(size_t(*&nodesTags)[constants::tetrahedron::N_FACES][constants::triangle::N_NODES], size_t& nNodesTags, const int entityTag);

            void getFaces(const int faceType, size_t* nodesTags, size_t nNodesTags, size_t** facesTags);
            void getTetrahedronsFaces(const size_t(*nodesTags)[constants::tetrahedron::N_FACES][constants::triangle::N_NODES], const size_t nNodesTags, size_t(*&facesTags)[constants::tetrahedron::N_FACES]);
            void getTetrahedronsFaces(size_t(*&nodesTags)[constants::tetrahedron::N_FACES][constants::triangle::N_NODES],
                                      size_t(*&facesTags)[constants::tetrahedron::N_FACES],
                                      size_t& nFaces,
                                      const int entityTag);

            void getTetrahedronsFacesOnSurface(const int surfaceTag, size_t*& facesTags, size_t& nFaces);



            void getNodesByElementType(const int elementType, size_t* &nodesTags, size_t &nNodesTags, double* &nodesCoordinates, const int tag = DEFAULT_TAG);

            void addElementsByType(const int tag, const int elementType, size_t* nodesTags, size_t nNodesTags, size_t* elementsTags = nullptr, size_t nElements = 0);
            void addTriangles(const int tag, size_t* nodesTags, size_t nNodesTags);

            void getNode(const size_t nodeTag, Coordinates& node);

            void get3Nodes(const size_t* nodeTagIt, Coordinates* nodeIt);
            
            void getNodes(const size_t* nodeTagIt, const size_t nNodes, Coordinates* nodeIt);

            void get3DElementsByCoordinates(const double x, const double y, const double z, size_t*& elementsTags, size_t& nElementsTags);
            
            int getElementEntityTag(const size_t elementTag);
            
            void getElement(const size_t elementTag, size_t* &nodesTags, size_t nNodes, int& tag);
            
        }
    }
};

