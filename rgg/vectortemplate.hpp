/*********************************************
Reactor Geometry Generator
Argonne National Laboratory

Vector template class definition.
*********************************************/
#ifndef __RGG_VECTORTEMPLATE_H__
#define __RGG_VECTORTEMPLATE_H__

#include <string>
#include <iostream>
#include <cassert>


// defines the vector template class
template <class T>
class CVector
#ifdef __RGGMETER
    : public CArrayBase
#endif
{
    private:                      
        int m_nRows;      // number of rows in the vector
        T   *m_pCells;    // address where the vector of 
                          // type T is stored
        void ErrorHandler (int,int nR=0) const;
                          // handles error conditions
        void Release ();  // similar to destructor
        std::string m_szName;       // vector name

    public:                     
        CVector ();                  // default constructor
        CVector (int);               // constructor
        CVector (const char*);       // constructor
        CVector (const char*, int);  // constructor
        CVector (const CVector<T>&); // copy constructor
        ~CVector ();                 // destructor
        void SetSize (int);     // sets the size of the vector 
                                // used with the default constructor

        // helper functions
        // gets the current size of the vector
        int GetSize () const;   
        // gets the vector name
        void GetName (std::string& szName) const; 

        // vector manipulations (mutator)
        void Set (T);           // sets the value of all
                                // elements of a vector
        void SetName (const std::string&); // sets the name of the vector

        // overloaded operators
        T& operator() (int);            // row access
        const T& operator() (int) const;// row access
        CVector<T>& operator= (const CVector<T>&);  // overloaded = operator
};

// =============== definitions ===========================================
template <class T>
CVector<T>::CVector ()
// ---------------------------------------------------------------------------
// Function: default ctor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_pCells = NULL; 
    m_nRows = 0;
}

template <class T>
CVector<T>::CVector (int n)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    size of the vector
// Output:   none
// ---------------------------------------------------------------------------
{
    m_pCells = NULL; 
    m_nRows = 0;
    SetSize (n);
}

template <class T>
CVector<T>::CVector (const char* szName)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector name
// Output:   none
// ---------------------------------------------------------------------------
{
    m_pCells = NULL; 
    m_nRows = 0;
	m_szName = szName;
}

template <class T>
CVector<T>::CVector (const char* szName, int n)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector name and size
// Output:   none
// ---------------------------------------------------------------------------
{
    m_pCells = NULL; 
    m_nRows = 0;
	m_szName = szName;
    SetSize (n);
}

template <class T>
CVector<T>::CVector (const CVector<T>& A)
// ---------------------------------------------------------------------------
// Function: copy ctor
// Input:    vector
// Output:   none
// ---------------------------------------------------------------------------
{
    m_pCells = NULL; 
    m_nRows = A.GetSize();
	m_szName = A.m_szName;
    SetSize (m_nRows);
    for (int i=1; i <= m_nRows; i++)
        m_pCells[i] = A.m_pCells[i];
}

template <class T>
void CVector<T>::SetSize (int nR)
// ---------------------------------------------------------------------------
// Function: dynamically allocates memory
// Input:    vector size
// Output:   none
// ---------------------------------------------------------------------------
{
    // check whether NR is legal
    if (nR <= 0) ErrorHandler (3);
    Release ();
    try {m_pCells = new T [nR + 1];}
    catch (std::bad_alloc) {ErrorHandler (1);}
    m_nRows = nR;
#ifdef __RGGMETER
    m_dAllocated += static_cast<double>(sizeof(T)*(nR+1));
#endif
}

template <class T>
CVector<T>::~CVector ()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // deallocate storage
    Release ();
}

template <class T>
void CVector<T>::Release ()
// ---------------------------------------------------------------------------
// Function: dynamically deallocates memory
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // deallocate storage
    if (m_pCells != NULL)
    {
        delete [] m_pCells;
#ifdef __RGGMETER
        m_dDeAllocated += static_cast<double>(sizeof(T)*(m_nRows+1));
#endif
        m_pCells = NULL; 
        m_nRows = 0;
    }
}

// =============== member functions ===========================================
template <class T>
void CVector<T>::SetName (const std::string& szName) 
// ---------------------------------------------------------------------------
// Function: sets the name of the vector
// Input:    vector name
// Output:   none
// ---------------------------------------------------------------------------
{
    m_szName = szName;
}

template <class T>
void CVector<T>::GetName (std::string& szName) const
// ---------------------------------------------------------------------------
// Function: gets the name of the vector
// Input:    string to hold vector name
// Output:   vector name
// ---------------------------------------------------------------------------
{
    szName = m_szName;
}

template <class T>
void CVector<T>::Set (T dV)
// ---------------------------------------------------------------------------
// Function: sets the value of all the elements in the vector to the
//           specified value
// Input:    specified value
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= m_nRows; i++) 
        m_pCells[i] = dV;
}

template <class T>
int CVector<T>::GetSize () const
// ---------------------------------------------------------------------------
// Function: gets the size of the vector
// Input:    none
// Output:   returns the size
// ---------------------------------------------------------------------------
{
    return m_nRows;
}

// ==================== Overloaded Operators ========================
#ifdef _DEBUG
template <class T>
T& CVector<T>::operator() (int nR)  // T& is reference
// ---------------------------------------------------------------------------
// Function: overloaded () operator to access vector contents
//           carries out bound checking
// Input:    index or location
// Output:   value at the specified index
// ---------------------------------------------------------------------------
{
    // row-column reference in bounds?
    if (nR <= 0 || nR > m_nRows)
    {
        ErrorHandler (2,nR);
        return m_pCells[1];
    }
    else
        return m_pCells[nR];
}
#else
template <class T>
inline T& CVector<T>::operator() (int nR)
// ---------------------------------------------------------------------------
// Function: overloaded () operator to access vector contents
// Input:    index or location
// Output:   value at the specified index
// ---------------------------------------------------------------------------
{
    return m_pCells[nR];
}
#endif

#ifdef _DEBUG
template <class T>
const T& CVector<T>::operator() (int nR) const
// ---------------------------------------------------------------------------
// Function: overloaded () operator to access vector contents
// Input:    index or location
// Output:   value at the specified index
// ---------------------------------------------------------------------------
{
    // row-column reference in bounds?
    if (nR <= 0 || nR > m_nRows) {
        ErrorHandler (2,nR);
        return m_pCells[1];
    }
    else
        return m_pCells[nR];
}
#else
template <class T>
inline const T& CVector<T>::operator() (int nR) const
// ---------------------------------------------------------------------------
// Function: overloaded () operator to access vector contents
// Input:    index or location
// Output:   value at the specified index
// ---------------------------------------------------------------------------
{
    return m_pCells[nR];
}
#endif

template <class T>
CVector<T>& CVector<T>::operator= (const CVector& matarg)
// ---------------------------------------------------------------------------
// Function: overloaded = operator 
// Input:    vector to use as rvalue
// Output:   modified values
// ---------------------------------------------------------------------------
{
    // check whether vector is assigned to itself
    if (this != &matarg)
    {
        // compatible vectors?
        if (m_nRows != matarg.m_nRows)
        {
            ErrorHandler (4);
            return *this;
        }
        // now copy
        for (int i=1; i <= matarg.m_nRows; i++)
            m_pCells[i] = matarg.m_pCells[i];
    }

    return *this;
}

// ==================== Error Handler ========================
template <class T>
void CVector<T>::ErrorHandler (int nErrorCode, int nR) const
// ---------------------------------------------------------------------------
// Function: channels error message via std:err
// Input:    error code and optional int value
// Output:   none
// ---------------------------------------------------------------------------
{
    std::cerr << "CVector::Vector Name: " << m_szName << ".\n";
    switch (nErrorCode)
    {
        case 1:
            std::cerr << "Memory allocation failure.\n";
        break;
        case 2:
            std::cerr << "Row index is out of bounds.\n";
        break;
        case 3:
            std::cerr << "Constructor. Invalid number of rows "
                         "or columns.\n";
        break;
        case 4:
            std::cerr << "Constructor. Incompatible vectors.\n";
        break;
    }
    // exit (1);
}

#endif
