#include "GMSHProxy.h"

#include <stdexcept>

extern "C" {
#include "gmshc.h"
}


namespace GMSHProxy
{
    inline static void throwLastError()
    {
        int ierr = 0;
        char* api_error_;
        gmshLoggerGetLastError(&api_error_, &ierr);
        if (ierr) throw std::runtime_error("Could not get last error");
        std::string error = std::string(api_error_);
        gmshFree(api_error_);
        throw std::runtime_error(error);
    }

    void initialize(int argc, char** argv,
                    const bool readConfigFiles,
                    const bool run)
    {
        int ierr = 0;
        gmshInitialize(argc, argv, (int)readConfigFiles, (int)run, &ierr);
        if (ierr)
        {
            throwLastError();
        }
    }

    void finalize()
    {
        int ierr;
        gmshFinalize(&ierr);
        if (ierr)
        {
            throwLastError();
        }
    }

    void open(const char modelName[])
    {
        int ierr;
        gmshOpen(modelName, &ierr);
        if (ierr)
        {
            throwLastError();
        }
    }

    void fltkRun()
    {
        int ierr;
        gmshFltkRun("", & ierr);
        if (ierr)
        {
            throwLastError();
        }
    }

    void free(void* p)
    {
        gmshFree(p);
    }

    namespace model
    {
        void extractTags(const int* dimTagIt, const size_t nTags, int* tagIt)
        {
            for (size_t i = 0; i < nTags; ++i)
            {
                ++dimTagIt;
                *tagIt = *dimTagIt;

                ++dimTagIt;
                ++tagIt;
            }
        }

        void getEntitiesForPhysicalName(const char name[], int* &dimTags, size_t &nDimTags)
        {
            int ierr = 0;
            gmshModelGetEntitiesForPhysicalName(name, &dimTags, &nDimTags, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void getEntitiesForPhysicalName(const char name[], int*& dimTags, unsigned int& nEntities)
        {
            int ierr = 0;

            size_t nDimTags;

            gmshModelGetEntitiesForPhysicalName(name, &dimTags, &nDimTags, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            nEntities = nDimTags / 2;
        }

        void getPhysicalGroups(int* &dimTags, size_t &nDimTags, int dim)
        {
            int ierr = 0;
            gmshModelGetPhysicalGroups(&dimTags, &nDimTags, dim, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void getPhysicalGroups(int*& dimTags, unsigned int& nGroups, int dim)
        {
            int ierr = 0;
            size_t nDimTags;
            gmshModelGetPhysicalGroups(&dimTags, &nDimTags, dim, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            nGroups = nDimTags / 2;
        }

        void getPhysicalGroups(const int dim, int** tags, size_t* nTags)
        {
            int* dimTags;
            size_t nDimTags;
            int ierr = 0;
            gmshModelGetPhysicalGroups(&dimTags, &nDimTags, dim, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            *nTags = nDimTags / 2;
            *tags = (int*)gmshMalloc(*nTags * sizeof(int));

            extractTags(dimTags, *nTags, *tags);
            gmshFree(dimTags);
        }

        void getPhysicalName(const int dim, const int tag, char* &physicalName)
        {
            int ierr = 0;
            gmshModelGetPhysicalName(dim, tag, &physicalName, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void getEntitiesForPhysicalGroup(const int dim, const int physicalTag, int* &entitiesTags, size_t &nEntities)
        {
            int ierr = 0;
            gmshModelGetEntitiesForPhysicalGroup(dim, physicalTag, &entitiesTags, &nEntities, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void getEntitiesForPhysicalGroup(const int dim, const int physicalTag, int*& entitiesTags, unsigned int& nEntities)
        {
            int ierr = 0;
            size_t nEntitiesLong;
            gmshModelGetEntitiesForPhysicalGroup(dim, physicalTag, &entitiesTags, &nEntitiesLong, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            nEntities = nEntitiesLong;
        }

        void getEntities(int* &dimTags, size_t &nDimTags, const int dim)
        {
            int ierr = 0;
            gmshModelGetEntities(&dimTags, &nDimTags, dim, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void getEntities(int*& dimTags, unsigned int& nEntities, const int dim)
        {
            int ierr = 0;
            size_t nDimTags;
            gmshModelGetEntities(&dimTags, &nDimTags, dim, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            nEntities = nDimTags / 2;
        }

        void getEntities(const int dim, int** tags, size_t* nEntities)
        {
            int* dimTags;
            size_t nDimTags;
            int ierr = 0;
            gmshModelGetEntities(&dimTags, &nDimTags, dim, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            *nEntities = nDimTags / 2;
            *tags = (int*)gmshMalloc(*nEntities * sizeof(int));

            extractTags(dimTags, *nEntities, *tags);
            gmshFree(dimTags);
        }

        void getSurfaceUpEntities(const int surfaceTag, int*& upEntitiesTags, size_t &nUpEntities)
        {
            int* downEntitiesTags;
            size_t nDownEntities;
            int ierr;
            gmshModelGetAdjacencies(constants::DIMENSION_2D, surfaceTag, &upEntitiesTags, &nUpEntities, &downEntitiesTags, &nDownEntities, &ierr);

            if (ierr)
            {
                throwLastError();
            }

            free(downEntitiesTags);
        }

        void getPhysicalGroupsForEntity(const int dim, const int entityTag, int* &physicalTags, size_t &nPhysicalGroups)
        {
            int ierr;
            gmshModelGetPhysicalGroupsForEntity(dim, entityTag, &physicalTags, &nPhysicalGroups, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void getBoundaries(const int entityDim, const int entityTag, int* &boundariesDimTags, size_t &nBoundariesDimTags)
        {
            int entityDimTag[] = { entityDim, entityTag };

            int ierr;
            gmshModelGetBoundary(entityDimTag, 2, &boundariesDimTags, &nBoundariesDimTags, FALSE_C, FALSE_C, FALSE_C, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void getBoundaries(const int entityDim, const int entityTag, int*& boundariesDimTags, unsigned int& nBoundaries)
        {
            int entityDimTag[] = { entityDim, entityTag };
            size_t nBoundariesDimTags;
            int ierr;
            gmshModelGetBoundary(entityDimTag, 2, &boundariesDimTags, &nBoundariesDimTags, FALSE_C, FALSE_C, FALSE_C, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            nBoundaries = nBoundariesDimTags / 2;
        }

        void getBoundariesFor3DEntity(const int entityTag, int*& boundariesTags, size_t& nBoundariesTags)
        {
            int entityDimTag[] = { constants::DIMENSION_3D, entityTag };

            int* boundariesDimTags;
            size_t nBoundariesDimTags;

            int ierr;
            gmshModelGetBoundary(entityDimTag, 2, &boundariesDimTags, &nBoundariesDimTags, FALSE_C, FALSE_C, FALSE_C, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            nBoundariesTags = nBoundariesDimTags / 2;

            boundariesTags = (int*)gmshMalloc(nBoundariesTags * sizeof(int));
            extractTags(boundariesDimTags, nBoundariesTags, boundariesTags);

            gmshFree(boundariesDimTags);
        }

        void getSurfaceBoundariesNodes(const int surfaceTag, int*& nodesTags, size_t& nNodes)
        {
            int surfaceDimTag[2] = { constants::DIMENSION_2D, surfaceTag };
            int* outDimTags;
            size_t nOutDimTags;

            int ierr;
            gmshModelGetBoundary(surfaceDimTag, 2, &outDimTags, &nOutDimTags, 0, 0, 1, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            size_t nTags = nOutDimTags / 2;
            nodesTags = (int*)gmshMalloc(nTags * sizeof(int));
            extractTags(outDimTags, nTags, nodesTags);

            free(outDimTags);
        }

        void getParametrizationOnSurface(const int surfaceTag, const double* coordinates, const size_t nCoordinates, double* &parametricCoordinates, size_t &nParametricCoordinates)
        {
            int ierr;
            gmshModelGetParametrization(constants::DIMENSION_2D, surfaceTag, coordinates, nCoordinates, &parametricCoordinates, &nParametricCoordinates, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void getNormalsByParametric(const int surfaceTag, const double* parametricCoordinates, const size_t nParametricCoords, double*& normalsCoordinates, size_t& nNormalsCoordinates)
        {
            int ierr;
            gmshModelGetNormal(surfaceTag, parametricCoordinates, nParametricCoords, &normalsCoordinates, &nNormalsCoordinates, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void getNormalsByCoordinates(const int surfaceTag, const double* coordinates, const size_t nCoordinates, double*& normalsCoordinates)
        {
            double* parametricCoordinates;
            size_t nParametricCoordinates;
            size_t nNormalsCoordinates;

            int ierr;
            gmshModelGetParametrization(constants::DIMENSION_2D, surfaceTag, coordinates, nCoordinates, &parametricCoordinates, &nParametricCoordinates, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            gmshModelGetNormal(surfaceTag, parametricCoordinates, nParametricCoordinates, &normalsCoordinates, &nNormalsCoordinates, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            free(parametricCoordinates);
        }

        void removeEntities(const int dim, const int* tagsIt, const size_t nTags)
        {
            int* dimTags = new int[2 * nTags];
            int* dimTagsIt = dimTags;
            for (size_t i = 0; i < nTags; ++i)
            {
                *dimTagsIt = dim;
                ++dimTagsIt;

                *dimTagsIt = *tagsIt;
                ++dimTagsIt;

                ++tagsIt;
            }

            int ierr;
            gmshModelRemoveEntities(dimTags, 2 * nTags, FALSE_C, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void getSurafceType(const int surfaceTag, char*& surfaceType)
        {
            int ierr;
            gmshModelGetType(constants::DIMENSION_2D, surfaceTag, &surfaceType, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        bool isSurfacePlane(const int surfaceTag)
        {
            char* surfaceType;

            int ierr;
            gmshModelGetType(constants::DIMENSION_2D, surfaceTag, &surfaceType, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            bool flag = strcmp(surfaceType, "Plane");
            gmshFree(surfaceType);

            return !flag;
        }

        int addDescreteEntity(int dim)
        {
            int ierr;
            int entityTag = gmshModelAddDiscreteEntity(dim, DEFAULT_TAG, nullptr, 0, &ierr);
            if (ierr)
            {
                throwLastError();
            }
            return entityTag;
        }

        int addDescreteEntity2D()
        {
            int ierr;
            int entityTag = gmshModelAddDiscreteEntity(constants::DIMENSION_2D, DEFAULT_TAG, nullptr, 0, &ierr);
            if (ierr)
            {
                throwLastError();
            }
            return entityTag;
        }

        void get2DPhysicalNames(const int* physicalTagsIt, const size_t nTags, char** physicalNamesIt)
        {
            int ierr = 0;
            for (size_t i = 0; i < nTags; ++i)
            {
                if (*physicalTagsIt == 0)
                {
                    *physicalNamesIt = nullptr;
                }
                else
                {
                    gmshModelGetPhysicalName(constants::DIMENSION_2D, *physicalTagsIt, physicalNamesIt, &ierr);
                    if (ierr)
                    {
                        throwLastError();
                    }
                }
                ++physicalTagsIt;
                ++physicalNamesIt;
            }
        }

        void deletePhysicalNames(char** physicalNamesIt, const size_t nPhysicalNames)
        {
            for (size_t i = 0; i < nPhysicalNames; ++i)
            {
                free(*physicalNamesIt);
                ++physicalNamesIt;
            }
            
        }

        void getCurrent(char*& modelName)
        {
            int ierr;
            gmshModelGetCurrent(&modelName, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void addNodes(const double* nodesCoordinateIt, const size_t nNodes, int* nodesTagsIt)
        {
            int ierr;
            for (size_t i = 0; i < nNodes; ++i)
            {
                *nodesTagsIt = gmshModelOccAddPoint(*nodesCoordinateIt, *(nodesCoordinateIt + 1), *(nodesCoordinateIt + 2), DEFAULT_MESH_SIZE, DEFAULT_TAG, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                ++nodesTagsIt;
                nodesCoordinateIt += 3;
            }
        }

        void setCurrent(const char modleName[])
        {
            int ierr;
            gmshModelSetCurrent(modleName , &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void remove()
        {
            int ierr;
            gmshModelRemove(&ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void removeEntities(const int* dimTags, const size_t nDimTags, bool recursive)
        {
            int ierr;
            gmshModelRemoveEntities(dimTags, nDimTags, recursive, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void remove2DEntities(const int* tags, const size_t nTags, bool recursive)
        {
            size_t nDimTags = 2 * nTags;
            int* dimTags = new int[nDimTags];
            int* dimTagsIt = dimTags;

            for (size_t i = 0; i < nTags; ++i)
            {
                *dimTagsIt = constants::DIMENSION_2D;
                ++dimTagsIt;
                *dimTagsIt = *tags;
                ++dimTagsIt;
            }

            int ierr;
            gmshModelRemoveEntities(dimTags, nDimTags, recursive, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            delete[] dimTags;
        }


        namespace mesh
        {
            void generate(int dim)
            {
                int ierr;
                gmshModelMeshGenerate(constants::DIMENSION_3D, &ierr);
                if (ierr)
                {
                   throwLastError();
                }
            }

            void createFaces(const int* dimTags, const size_t nDimTags)
            {
                int ierr;
                gmshModelMeshCreateFaces(dimTags, nDimTags, &ierr);
                if (ierr)
                {
                    throwLastError();
                }
            }

            void createFaces()
            {
                int ierr;
                gmshModelMeshCreateFaces(nullptr, 0, &ierr);
                if (ierr)
                {
                    throwLastError();
                }
            }

            void getJacobians(const int elementType, double** jacobians, double** determinants, double** coord, size_t* nElements, const int tag)
            {
                const double localCoord[LocalCoordinates3D::COUNT] = { 0.0, 0.0, 0.0 };
                int ierr;
                size_t jacobians_n, coord_n;
                gmshModelMeshGetJacobians(elementType,
                                          localCoord,
                                          LocalCoordinates3D::COUNT,
                                          jacobians,
                                          &jacobians_n,
                                          determinants,
                                          nElements,
                                          coord,
                                          &coord_n,
                                          tag,
                                          TASK_DEFAULT,
                                          N_TASK_DEFAULT,
                                          &ierr);

                if (ierr)
                {
                    throwLastError();
                }
            }

            void getJacobians(const int elementType, double* &jacobians, double* &determinants, double* &coord, size_t &nElements, const int tag)
            {
                const double localCoord[LocalCoordinates3D::COUNT] = { 0.0, 0.0, 0.0 };
                int ierr;
                size_t jacobians_n, coord_n;
                gmshModelMeshGetJacobians(elementType,
                    localCoord,
                    LocalCoordinates3D::COUNT,
                    &jacobians,
                    &jacobians_n,
                    &determinants,
                    &nElements,
                    &coord,
                    &coord_n,
                    tag,
                    TASK_DEFAULT,
                    N_TASK_DEFAULT,
                    &ierr);

                if (ierr)
                {
                    throwLastError();
                }
            }

            void getJacobians(const int elementType, double (*&transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT], double* &determinants, Coordinates* &initPoints, size_t &nElements, const int tag)
            {
                const double localCoord[LocalCoordinates3D::COUNT] = { 0.0, 0.0, 0.0 };
                int ierr;
                double* transpJacobianMatrixesElements;
                double* initPointsCoordinates;
                size_t jacobians_n, coord_n;
                gmshModelMeshGetJacobians(elementType,
                    localCoord,
                    LocalCoordinates3D::COUNT,
                    &transpJacobianMatrixesElements,
                    &jacobians_n,
                    &determinants,
                    &nElements,
                    &initPointsCoordinates,
                    &coord_n,
                    tag,
                    TASK_DEFAULT,
                    N_TASK_DEFAULT,
                    &ierr);

                transpJacobianMatrixes = (double(*)[Coordinates::COUNT * LocalCoordinates3D::COUNT])transpJacobianMatrixesElements;
                initPoints = (Coordinates*)initPointsCoordinates;

                if (ierr)
                {
                    throwLastError();
                }
            }

            void getInitJacobian(const int elementTag, double*& jacobians, double &determinant)
            {
                int ierr;
                double localCoordinates[LocalCoordinates3D::COUNT] = { 0.0, 0.0, 0.0 };
                double* coordinates, *determinants;
                size_t nCoordinates, nDeterminants, nJacobianElements;
                gmshModelMeshGetJacobian(elementTag, localCoordinates, LocalCoordinates3D::COUNT, &jacobians, &nJacobianElements, &determinants, &nDeterminants, &coordinates, &nCoordinates, &ierr);

                if (ierr)
                {
                    throwLastError();
                }

                determinant = *determinants;

                free(coordinates);
                free(determinants);
            }

            void getTriangleJacobians(double(*&transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                                      double*& determinants,
                                      Coordinates*& initPoints,
                                      const int tag)
            {
                const double localCoord[LocalCoordinates3D::COUNT] = { 0.0, 0.0, 0.0 };
                int ierr;
                double* transpJacobianMatrixesElements;
                double* initPointsCoordinates;
                size_t jacobians_n, coord_n, nElements;
                gmshModelMeshGetJacobians(TRIANGLE_TYPE,
                                          localCoord,
                                          LocalCoordinates3D::COUNT,
                                          &transpJacobianMatrixesElements,
                                          &jacobians_n,
                                          &determinants,
                                          &nElements,
                                          &initPointsCoordinates,
                                          &coord_n,
                                          tag,
                                          TASK_DEFAULT,
                                          N_TASK_DEFAULT,
                                          &ierr);

                transpJacobianMatrixes = (double(*)[Coordinates::COUNT * LocalCoordinates3D::COUNT])transpJacobianMatrixesElements;
                initPoints = (Coordinates*)initPointsCoordinates;

                if (ierr)
                {
                    throwLastError();
                }
            }

            void getTetrahedronsJacobians(double(*&transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                                         double*& determinants,
                                         Coordinates*& initPoints,
                                         size_t& nElements,
                                         const int tag)
            {
                const double localCoord[LocalCoordinates3D::COUNT] = { 0.0, 0.0, 0.0 };
                int ierr;
                double* transpJacobianMatrixesElements;
                double* initPointsCoordinates;
                size_t jacobians_n, coord_n;
                gmshModelMeshGetJacobians(TETRAHEDRON_TYPE,
                                          localCoord,
                                          LocalCoordinates3D::COUNT,
                                          &transpJacobianMatrixesElements,
                                          &jacobians_n,
                                          &determinants,
                                          &nElements,
                                          &initPointsCoordinates,
                                          &coord_n,
                                          tag,
                                          TASK_DEFAULT,
                                          N_TASK_DEFAULT,
                                          &ierr);

                transpJacobianMatrixes = (double(*)[Coordinates::COUNT * LocalCoordinates3D::COUNT])transpJacobianMatrixesElements;
                initPoints = (Coordinates*)initPointsCoordinates;

                if (ierr)
                {
                    throwLastError();
                }
            }

            void getTetrahedronsJacobians(double(*&transpJacobianMatrixes)[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                                         double*& determinants,
                                         Coordinates*& initPoints,
                                         const int tag)
            {
                int ierr;

                const double localCoord[LocalCoordinates3D::COUNT] = { 0.0, 0.0, 0.0 };
                double* transpJacobianMatrixesElements;
                double* initPointsCoordinates;
                size_t nJacobians, nDeterminants, nCoordinates;
                gmshModelMeshGetJacobians(TETRAHEDRON_TYPE,
                                          localCoord,
                                          LocalCoordinates3D::COUNT,
                                          &transpJacobianMatrixesElements,
                                          &nJacobians,
                                          &determinants,
                                          &nDeterminants,
                                          &initPointsCoordinates,
                                          &nCoordinates,
                                          tag,
                                          TASK_DEFAULT,
                                          N_TASK_DEFAULT,
                                          &ierr);

                transpJacobianMatrixes = (double(*)[Coordinates::COUNT * LocalCoordinates3D::COUNT])transpJacobianMatrixesElements;
                initPoints = (Coordinates*)initPointsCoordinates;

                if (ierr)
                {
                    throwLastError();
                }
            }

            void getTrianglesNodes(size_t(*&traianglesNodesTags)[constants::triangle::N_NODES], Coordinates(*&traianglesNodes)[constants::triangle::N_NODES], const int tag)
            {
                size_t* tags;
                size_t nTags;

                double* coordinates;
                size_t nCoordinates;

                double* parametricCoordinates;
                size_t nParametricCoordinates;

                int ierr;

                gmshModelMeshGetNodesByElementType(TRIANGLE_TYPE, &tags, &nTags, &coordinates, &nCoordinates, &parametricCoordinates, &nParametricCoordinates, tag, FALSE_C, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                traianglesNodesTags = (size_t(*)[constants::triangle::N_NODES])tags;
                traianglesNodes = (Coordinates(*)[constants::triangle::N_NODES])coordinates;
            }

            void getIntegrationPoints(const int elementType, const char integrationType[], double** localCoord, double** weights, uint8_t* nIntegrationSteps)
            {
                int ierr;
                size_t nLocalCoords;
                size_t nIntegrationStepsProxy;
                gmshModelMeshGetIntegrationPoints(elementType, integrationType, localCoord, &nLocalCoords, weights, &nIntegrationStepsProxy, &ierr);

                *nIntegrationSteps = nIntegrationStepsProxy;

                if (ierr)
                {
                    throwLastError();
                }
            }

            void getElementsByType(const int elementType, size_t* &elementTags, size_t &nElement, size_t* &nodeTags, const int tag, const size_t task, const size_t nTasks)
            {
                int ierr;
                size_t nNodes;
                gmshModelMeshGetElementsByType(elementType, &elementTags, &nElement, &nodeTags, &nNodes, tag, task, nTasks, &ierr);

                if (ierr)
                {
                    throwLastError();
                }
            }
            void getTetrahedrons(size_t*& tetrahedronsTags, size_t& nTetrahedrons, size_t(*&tetrahedronsNodesTags)[constants::tetrahedron::N_NODES], const int tag, const size_t task, const size_t nTasks)
            {
                int ierr;
                size_t nNodes;
                size_t* nodesTags;
                gmshModelMeshGetElementsByType(TETRAHEDRON_TYPE, &tetrahedronsTags, &nTetrahedrons, &nodesTags, &nNodes, tag, task, nTasks, &ierr);

                if (ierr)
                {
                    throwLastError();
                }

                tetrahedronsNodesTags = (size_t(*)[constants::tetrahedron::N_NODES])nodesTags;
            }


            void getFacesByElements(const int elementType, const int faceType, size_t** nodesTags, size_t* nNodesTags, const int entityTag)
            {
                int ierr;
                gmshModelMeshGetElementFaceNodes(elementType, faceType, nodesTags, nNodesTags, entityTag, FALSE_C, TASK_DEFAULT, N_TASK_DEFAULT, &ierr);
                if (ierr)
                {
                    throwLastError();
                }
            }

            void getTetrahedronsFacesByElements(size_t (*&nodesTags)[constants::tetrahedron::N_FACES][constants::triangle::N_NODES], size_t &nNodesTags, const int entityTag)
            {
                int ierr;
                gmshModelMeshGetElementFaceNodes(TETRAHEDRON_TYPE, TRIANGLE_FACE_TYPE, (size_t**)(&nodesTags), &nNodesTags, entityTag, FALSE_C, TASK_DEFAULT, N_TASK_DEFAULT, &ierr);
                if (ierr)
                {
                    throwLastError();
                }
            }

            void getFaces(const int faceType, size_t* nodesTags, size_t nNodesTags, size_t** facesTags)
            {
                int ierr;
                int* facesOrientations;
                size_t nFacesOrientations, nFaces;

                gmshModelMeshGetFaces(faceType, nodesTags, nNodesTags, facesTags, &nFaces, &facesOrientations, &nFacesOrientations, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                gmshFree(facesOrientations);
            }

            void getTetrahedronsFaces(const size_t(*nodesTags)[constants::tetrahedron::N_FACES][constants::triangle::N_NODES], const size_t nNodesTags, size_t (*&facesTags)[constants::tetrahedron::N_FACES])
            {
                int ierr;
                int* facesOrientations;
                size_t nFacesOrientations, nFaces;

                gmshModelMeshGetFaces(TRIANGLE_FACE_TYPE, (const size_t*)nodesTags, nNodesTags, (size_t**)(&facesTags), &nFaces, &facesOrientations, &nFacesOrientations, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                gmshFree(facesOrientations);
            }

            void getTetrahedronsFaces(size_t(*&nodesTags)[constants::tetrahedron::N_FACES][constants::triangle::N_NODES],
                                      size_t(*&facesTags)[constants::tetrahedron::N_FACES],
                                      size_t &nFaces,
                                      const int entityTag)
            {
                int ierr;
                size_t nNodes;
                gmshModelMeshGetElementFaceNodes(TETRAHEDRON_TYPE, TRIANGLE_FACE_TYPE, (size_t**)(&nodesTags), &nNodes, entityTag, FALSE_C, TASK_DEFAULT, N_TASK_DEFAULT, &ierr);

                int* facesOrientations;
                size_t nFacesOrientations;

                gmshModelMeshGetFaces(TRIANGLE_FACE_TYPE, (const size_t*)nodesTags, nNodes, (size_t**)(&facesTags), &nFaces, &facesOrientations, &nFacesOrientations, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                gmshFree(facesOrientations);
            }

            void getTetrahedronsFacesOnSurface(const int surfaceTag, size_t* &facesTags, size_t& nFaces)
            {
                int ierr;

                size_t* elementTags;
                size_t* nodeTags;
                size_t nElements;
                size_t nNodes;

                gmshModelMeshGetElementsByType(TRIANGLE_TYPE, &elementTags, &nElements, &nodeTags, &nNodes, surfaceTag, TASK_DEFAULT, N_TASK_DEFAULT, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                gmshFree(elementTags);

                int* facesOrientations;
                size_t nFacesOrientations;
                gmshModelMeshGetFaces(TRIANGLE_FACE_TYPE, (const size_t*)nodeTags, nNodes, (size_t**)(&facesTags), &nFaces, &facesOrientations, &nFacesOrientations, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                gmshFree(facesOrientations);
            }

            void getNodesByElementType(const int elementType, size_t* &nodesTags, size_t &nNodesTags, double* &nodesCoordinates, const int tag)
            {
                size_t nCoordinates;
                double* parametricCoordinates;
                size_t nParametricCoordinates;
                int ierr;
                gmshModelMeshGetNodesByElementType(elementType, &nodesTags, &nNodesTags, &nodesCoordinates, &nCoordinates, &parametricCoordinates, &nParametricCoordinates, tag, FALSE_C, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

            }

            void addElementsByType(const int tag, const int elementType, size_t* nodesTags, size_t nNodesTags, size_t* elementsTags, size_t nElements)
            {
                int ierr;
                gmshModelMeshAddElementsByType(tag, elementType, elementsTags, nElements, nodesTags, nNodesTags, &ierr);
                if (ierr)
                {
                    throwLastError();
                }
            }

            void addTriangles(const int tag, size_t* nodesTags, size_t nNodesTags)
            {
                int ierr;
                gmshModelMeshAddElementsByType(tag, TRIANGLE_TYPE, nullptr, 0, nodesTags, nNodesTags, &ierr);
                if (ierr)
                {
                    throwLastError();
                }
            }

            void getNode(const size_t nodeTag, Coordinates& node)
            {
                double* nodeCoordinates;
                size_t nCoordinates;
                double* parametricCoordinates;
                size_t nParametricCoordinates;

                int dim, tag, ierr;

                gmshModelMeshGetNode(nodeTag, &nodeCoordinates, &nCoordinates, &parametricCoordinates, &nParametricCoordinates, &dim, &tag, &ierr);

                if (ierr)
                {
                    throwLastError();
                }

                node.x = nodeCoordinates[0];
                node.y = nodeCoordinates[1];
                node.z = nodeCoordinates[2];

                gmshFree(nodeCoordinates);
                gmshFree(parametricCoordinates);
            }

            void get3Nodes(const size_t* nodeTagIt, Coordinates* nodeIt)
            {
                double* nodeCoordinates;
                size_t nCoordinates;
                double* parametricCoordinates;
                size_t nParametricCoordinates;

                int dim, tag, ierr;

                for (size_t i = 0; i < 3; ++i)
                {
                    gmshModelMeshGetNode(*nodeTagIt, &nodeCoordinates, &nCoordinates, &parametricCoordinates, &nParametricCoordinates, &dim, &tag, &ierr);

                    if (ierr)
                    {
                        throwLastError();
                    }


                    nodeIt->x = nodeCoordinates[0];
                    nodeIt->y = nodeCoordinates[1];
                    nodeIt->z = nodeCoordinates[2];

                    gmshFree(nodeCoordinates);
                    gmshFree(parametricCoordinates);

                    ++nodeTagIt;
                    ++nodeIt;
                }
            }

            void getNodes(const size_t* nodeTagIt, const size_t nNodes, Coordinates* nodeIt)
            {
                double* nodeCoordinates;
                size_t nCoordinates;
                double* parametricCoordinates;
                size_t nParametricCoordinates;

                int dim, tag, ierr;

                for (size_t i = 0; i < nNodes; ++i)
                {
                    gmshModelMeshGetNode(*nodeTagIt, &nodeCoordinates, &nCoordinates, &parametricCoordinates, &nParametricCoordinates, &dim, &tag, &ierr);

                    if (ierr)
                    {
                        throwLastError();
                    }


                    nodeIt->x = nodeCoordinates[0];
                    nodeIt->y = nodeCoordinates[1];
                    nodeIt->z = nodeCoordinates[2];

                    gmshFree(nodeCoordinates);
                    gmshFree(parametricCoordinates);

                    ++nodeTagIt;
                    ++nodeIt;
                }
            }

            void extractTags(const int* dimTagsIt, const size_t nTags, int* tagIt)
            {
                for (size_t i = 0; i < nTags; ++i)
                {
                    ++dimTagsIt;
                    *tagIt = *dimTagsIt;

                    ++dimTagsIt;
                    ++tagIt;
                }
            }

            void get3DElementsByCoordinates(const double x, const double y, const double z, size_t* &elementsTags, size_t& nElementsTags)
            {
                int ierr;
                gmshModelMeshGetElementsByCoordinates(x, y, z, &elementsTags, &nElementsTags, constants::DIMENSION_3D, FALSE_C, &ierr);
                if (ierr)
                {
                    throwLastError();
                }
            }

            int getElementEntityTag(const size_t elementTag)
            {
                int entityDim, entityTag;
                size_t* nodesTags;
                size_t nNodes;
                int elementType;

                int ierr;
                gmshModelMeshGetElement(elementTag, &elementType, &nodesTags, &nNodes, &entityDim, &entityTag, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                return entityTag;
            }

            void getElement(const size_t elementTag, size_t*& nodesTags, size_t nNodes, int& entityTag)
            {
                int elementType, entityDim, ierr;
                gmshModelMeshGetElement(elementTag, &elementType, &nodesTags, &nNodes, &entityDim, &entityTag, &ierr);
                if (ierr)
                {
                    throwLastError();
                }
            }
            
        }
    }
}; 

