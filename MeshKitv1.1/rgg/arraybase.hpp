/*********************************************
Dec,09
Reactor Geometry Generator
Argonne National Laboratory

Utility Library Function
*********************************************/
#ifndef __RGG_ARRAYBASE_H__
#define __RGG_ARRAYBASE_H__

class CArrayBase
{
    public:
        CArrayBase ();
        ~CArrayBase ();

        void static ShowStatistics ();

    protected:
        static double m_dAllocated;
        static double m_dDeAllocated;
};

#endif
