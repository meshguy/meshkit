/*
 * AF2AlgorithmResult.hpp
 *
 * \brief An AF2AlgorithmResult stores the result of executing
 * a two-dimensional advancing front algorithm.
 *
 * The result may be successful or unsuccessful.
 */
#ifndef AF2ALGORITHMRESULT_HPP
#define AF2ALGORITHMRESULT_HPP

// C++
#include <list>

// MeshKit
#include "meshkit/AF2Point3D.hpp"
#include "meshkit/AF2Polygon3D.hpp"

class AF2AlgorithmResult
{
  private:

    bool succeeded;
    const std::list<AF2Point3D*> points;
    const std::list<const AF2Polygon3D*> faces;

  public:

    /**
     * \brief Constructor for an unsuccessful result
     *
     * A constructor for an algorithm result that indicates the algorithm
     * terminated without successfully completing.
     */
    AF2AlgorithmResult();

    /**
     * \brief Constructor for a successful result
     *
     * A constructor to contain the results of a successfully completed
     * two-dimensional advancing front algorithm.
     *
     * This class takes ownership of the AF2Point3D and AF2Polygon3D
     * objects that are in the lists passed to the constructor.  The
     * aforementioned points and faces are deleted when this class
     * is destructed.
     *
     * \param pointsArg a list of the points that are in the mesh at
     *   the successful conclusion of the advancing front algorithm
     * \param facesArg a list of the faces that are in the mesh at
     *   the successful conclusion of the advancing front algorithm
     */
    AF2AlgorithmResult(const std::list<AF2Point3D*> & pointsArg,
        const std::list<const AF2Polygon3D*> & facesArg);

    /**
     * \brief Destructor
     */
    ~AF2AlgorithmResult();

    /**
     * \brief Copy constructor (throws exception)
     *
     * This object does not currently support copying, so this method
     * is implemented to throw an exception.
     */
    AF2AlgorithmResult(const AF2AlgorithmResult & toCopy);

    /**
     * \brief Assignment operator (throws exception)
     *
     * This object does not currently support assignment, so this method
     * is implemented to throw an exception.
     */
    AF2AlgorithmResult& operator=(const AF2AlgorithmResult & rhs);

    /**
     * \brief Get a pointer to the list of faces that are
     *   in the advancing front algorithm result.
     *
     * If the result was unsuccessful, this method will return NULL.
     *
     * \return the list of faces added to the mesh, if the result was
     *   successful
     */
    const std::list<const AF2Polygon3D*>* getFaces() const;

    /**
     * \brief Get a pointer to the list of points that are
     *   in the advancing front algorithm result.
     *
     * If the result was unsuccessful, this method will return NULL.
     *
     * \return the list of points added to the mesh, if the result was
     *   successful
     */
    const std::list<AF2Point3D*>* getPoints() const;

    /**
     * \brief Report whether the algorithm succeeded
     *
     * \return true if the algorithm succeeded, false otherwise
     */
    bool isSuccessful() const;
};

#endif
