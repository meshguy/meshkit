/*
 * AF2Neighborhood.hpp
 *
 * TODO: Document
 */
#ifndef AF2NEIGHBORHOOD_HPP
#define AF2NEIGHBORHOOD_HPP

// C++
#include <list>

// MeshKit
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2Edge2D.hpp"

class AF2Neighborhood
{

  public:

    /**
     * \brief TODO: Document method briefly
     *
     * TODO: Document method
     *
     * \param TODO: Document method parameters
     */
    const std::list<AF2Edge2D*>* getEdges2D() const;

    /**
     * \brief TODO: Document method briefly
     *
     * TODO: Document method
     *
     * \param TODO: Document method parameters
     */
    const std::list<AF2Point2D*>* getPoints2D() const;
};

#endif
