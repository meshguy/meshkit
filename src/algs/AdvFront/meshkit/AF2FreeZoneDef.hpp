/*
 * AF2FreeZoneDef.hpp
 *
 * A free zone definition is an object that can be used to define an
 * AF2FreeZone based on a binding of rule reference vertices to actual
 * coordinates and a level of quality, also known as a quality class,
 * at which the free zone is being applied.
 */

#ifndef AF2FREEZONEDEF_HPP
#define AF2FREEZONEDEF_HPP

// MeshKit
#include "meshkit/AF2FreeZone.hpp"
#include "meshkit/AF2Binding.hpp"

class AF2FreeZoneDef
{
  public:

    /**
     * \brief Make and return an independent copy of yourself.
     *
     * This method supports making copies of the concrete derived class that
     * an AF2FreeZoneDef pointer references.
     *
     * Deletion of the object pointed to by the pointer returned from this
     * method must be managed by the context that calls the method.
     */
    virtual AF2FreeZoneDef* clone() const = 0;

    /**
     * \brief Create an AF2FreeZone based on a binding and a quality class.
     *
     * This method allocates the AF2FreeZone and returns a pointer to it.
     * It is the responsibility of the calling context to make sure that
     * the AF2FreeZone is deleted when it is no longer needed.
     *
     * \param vertexBinding a binding of a rule's existing vertices to
     *   actual points that will be used to build a free zone relative
     *   to the actual points
     * \param qualityClass an integer strictly greater than zero
     *   designating the quality class of the free zone this method
     *   will make.  Lower numbers are higher quality and may produce
     *   larger free zones, depending on the implementation.
     */
    virtual AF2FreeZone* makeFreeZone(
        AF2Binding const & vertexBinding, int qualityClass) const = 0;
};

#endif
