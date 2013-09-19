/*!
 * \class matrixtemplate
 * \brief for engineering style matrix starting from 1
 *
 **/
#ifndef __MESHKIT_MATRIXTEMPLATE_H__
#define __MESHKIT_MATRIXTEMPLATE_H__

#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <stdlib.h>

// defines the Matrix class
template <class T>
class CMatrix
#ifdef __RGGMETER
    : public CArrayBase
#endif
{
    private:                      
        T   **m_pCells;   // address where the matrix of 
                          // type T is stored
        int m_nRows;      // number of rows in the matrix
        int m_nColumns;   // number of columns in the matrix
        std::string m_szName; // matrix name
        void ErrorHandler (int,int nR=0, int nC=0) const;
                                // handles error conditions
        void Release ();        // similar to destructor

    public:                     
        CMatrix ();                       // default constructor
        CMatrix (int, int);               // constructor
        CMatrix (const char *);           // constructor
        CMatrix (const char *, int, int); // constructor
        CMatrix (const CMatrix<T>&); // copy constructor
        ~CMatrix ();                 // destructor
        void SetSize (int, int);     // sets the size of the matrix
                                     // used with the default constructor

        // helper functions
        int GetRows () const;       // gets the current number of rows
        int GetColumns () const;    // gets the current number of columns
        void GetName (std::string&) const; // gets the matrix name

        // matrix manipulations (mutator)
        void Set (T);               // sets the value of all elements
                                    // of a matrix
        void SetName (const std::string&); // sets the matrix name
        T& operator() (int, int);               // row-col access
        const T& operator() (int, int) const;   // row-col access
        T& operator= (const CMatrix&);  // overloaded = operator
        int Read (std::ifstream& IFile);
};

// =============== definitions ===========================================
template <class T>
CMatrix<T>::CMatrix ()
// ---------------------------------------------------------------------------
// Function: default ctor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_pCells = NULL; 
    m_nRows = 0;
    m_nColumns = 0;
}

template <class T>
CMatrix<T>::CMatrix (int nR, int nC)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    # of rows and columns
// Output:   none
// ---------------------------------------------------------------------------
{
    m_pCells = NULL; 
    SetSize (nR, nC);
}

template <class T>
CMatrix<T>::CMatrix (const char* szName)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    matrix name
// Output:   none
// ---------------------------------------------------------------------------
{
    m_pCells = NULL; 
	m_szName = szName;
    m_nRows = 0;
    m_nColumns = 0;
}

template <class T>
CMatrix<T>::CMatrix (const char* szName, int nR, int nC)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    matrix name, # of rows and columns
// Output:   none
// ---------------------------------------------------------------------------
{
    m_pCells = NULL; 
	m_szName = szName;
    SetSize (nR, nC);
}

template <class T>
CMatrix<T>::CMatrix (const CMatrix<T>& A)
// ---------------------------------------------------------------------------
// Function: copy ctor
// Input:    matrix
// Output:   none
// ---------------------------------------------------------------------------
{
    m_pCells = NULL;
    m_nRows = A.m_nRows;
    m_nColumns = A.m_nColumns;
	m_szName = A.m_szName;
    SetSize (m_nRows, m_nColumns);
    for (int i=1; i <= m_nRows; i++)
    {
        for (int j=1; j <= m_nColumns; j++)
        {
            m_pCells[i][j] = A.m_pCells[i][j];
        }
    }
}

template <class T>
void CMatrix<T>::SetSize (int nR, int nC)
// ---------------------------------------------------------------------------
// Function: dynamically allocates memory
// Input:    matrix size (# of rows and columns)
// Output:   none
// ---------------------------------------------------------------------------
{
    // check whether nR and nC are legal
    if (nR <= 0 || nC <= 0)
        ErrorHandler (3);
    Release ();
    int size = nR*nC + 1;
    try {m_pCells = new T *[nR + 1];}
    catch (std::bad_alloc) {ErrorHandler (1);}
    try {m_pCells[0] = new T[size];}
    catch (std::bad_alloc) {ErrorHandler (1);}
    m_pCells[1] = m_pCells[0];
    for (int i=2; i <= nR; i++)
         m_pCells[i] = m_pCells[i-1]+nC;
    m_nRows = nR;
    m_nColumns = nC;
#ifdef __RGGMETER
    m_dAllocated += static_cast<double>(sizeof(T*)*(nR+1));
    m_dAllocated += static_cast<double>(sizeof(T)*size);
#endif
}

template <class T>
CMatrix<T>::~CMatrix ()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // deallocate storage
    if (m_pCells != NULL)
    {
        delete [] m_pCells[0];
        m_pCells[0] = NULL;
#ifdef __RGGMETER
        m_dDeAllocated += static_cast<double>(sizeof(T*)*(m_nRows+1));
#endif
    }
    if (m_pCells != NULL)
    {
        delete [] m_pCells;
        m_pCells = NULL;
#ifdef __RGGMETER
        int nSize = m_nRows*m_nColumns+1;
        m_dDeAllocated += static_cast<double>(sizeof(T)*nSize);
#endif
    }
    m_nRows = 0;
    m_nColumns = 0;
}

template <class T>
void CMatrix<T>::Release ()
// ---------------------------------------------------------------------------
// Function: dynamically deallocates memory
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // deallocate storage
    if (m_pCells != NULL)
    {
        delete [] m_pCells[0];
#ifdef __RGGMETER
        m_dDeAllocated += static_cast<double>(sizeof(T*)*(m_nRows+1));
#endif
    }
    if (m_pCells != NULL)
    {
        delete [] m_pCells;
#ifdef __RGGMETER
        int nSize = m_nRows*m_nColumns+1;
        m_dDeAllocated += static_cast<double>(sizeof(T)*nSize);
#endif
    }
    m_pCells = NULL; 
    m_nRows = 0;
    m_nColumns = 0;
}

// =============== member functions ===========================================
template <class T>
void CMatrix<T>::SetName (const std::string& szName) 
// ---------------------------------------------------------------------------
// Function: sets the name of the matrix
// Input:    matrix name
// Output:   none
// ---------------------------------------------------------------------------
{
    m_szName = szName;
}

template <class T>
void CMatrix<T>::GetName (std::string& szName) const
// ---------------------------------------------------------------------------
// Function: gets the name of the matrix
// Input:    string to hold matrix name
// Output:   matrix name
// ---------------------------------------------------------------------------
{
    szName = m_szName;
}

template <class T>
int CMatrix<T>::GetRows () const
// ---------------------------------------------------------------------------
// Function: gets the # of rows in the matrix
// Input:    none
// Output:   # of rows
// ---------------------------------------------------------------------------
{
    return (m_nRows);
}

template <class T>
int CMatrix<T>::GetColumns () const
// ---------------------------------------------------------------------------
// Function: gets the # of columns in the matrix
// Input:    none
// Output:   # of columns
// ---------------------------------------------------------------------------
{
    return (m_nColumns);
}

template <class T>
void CMatrix<T>::Set (T dV)
// ---------------------------------------------------------------------------
// Function: sets the value of all the elements in the matrix to the
//           specified value
// Input:    specified value
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= m_nRows; i++)
    {
        for (int j=1; j <= m_nColumns; j++)
        {
            m_pCells[i][j] = dV; // or, (*this)(i,j) = dV;
        }
    }
}

// ==================== Overloaded Operators ========================
template <class T>
T& CMatrix<T>::operator() (int nR, int nC)
// ---------------------------------------------------------------------------
// Function: overloaded () operator to access matrix contents
//           carries out bound checking
// Input:    row and column indices
// Output:   value at the specified indices
// ---------------------------------------------------------------------------
{
#ifdef _DEBUG
    if (nR <= 0 || nR > m_nRows || nC <= 0 || nC > m_nColumns)
    {
        ErrorHandler (2,nR,nC);
        return m_pCells[1][1];
    }
    else
        return m_pCells[nR][nC];
#else
        return m_pCells[nR][nC];
#endif
}

template <class T>
const T& CMatrix<T>::operator() (int nR, int nC) const
// ---------------------------------------------------------------------------
// Function: overloaded () operator to access matrix contents
//           carries out bound checking
// Input:    row and column indices
// Output:   value at the specified indices
// ---------------------------------------------------------------------------
{
#ifdef _DEBUG
    if (nR <= 0 || nR > m_nRows || nC <= 0 || nC > m_nColumns)
    {
        ErrorHandler (2,nR,nC);
        return m_pCells[1][1];
    }
    else
        return m_pCells[nR][nC];
#else
        return m_pCells[nR][nC];
#endif
}

template <class T>
T& CMatrix<T>::operator= (const CMatrix& matarg)
// ---------------------------------------------------------------------------
// Function: overloaded = operator 
// Input:    matrix to use as rvalue
// Output:   modified values
// ---------------------------------------------------------------------------
{
    // check whether matrix is assigned to itself
    if (this != &matarg)
    {
        // compatible matrices?
        if (m_nRows != matarg.m_nRows || m_nColumns != matarg.m_nColumns)
        {
            ErrorHandler (4);
            return (T&)(*this);
        }
        // now copy
        for (int i=1; i <= matarg.m_nRows; i++)
        {
            for (int j=1; j <= matarg.m_nColumns; j++)
            {
                m_pCells[i][j]= matarg.m_pCells[i][j];
            }
        }
    }

    return (T&)(*this);
}

template <class T>
int CMatrix<T>::Read (std::ifstream& IFile)
// ---------------------------------------------------------------------------
// Function: Reads a matrix rowwise
//    Input: input stream to read from
//   Output: the modified matrix
// ---------------------------------------------------------------------------
{
	// read matrix size
	int m_nRowsI, nColsI;
	IFile >> m_nRowsI >> nColsI;
	if (m_nRowsI <= 0 || nColsI <= 0)
		return 1;
	if (IFile.eof() || IFile.fail())
		return 1;

	// reallocate?
	if (m_nRows != m_nRowsI || m_nColumns != nColsI)
		SetSize (m_nRowsI, nColsI);
	
	T v;

	for (int i=1; i <= m_nRows; i++)
	{
		for (int j=1; j <= m_nColumns; j++)
		{
			IFile >> v;
			m_pCells[i][j] = v;
			if (IFile.eof() || IFile.fail())
				return 1;
		}
	}

	return 0;
}

// ==================== Error Handler ========================
template <class T>
void CMatrix<T>::ErrorHandler (int nErrorCode, int nR, int nC) const
// ---------------------------------------------------------------------------
// Function: channels error message via std:err
// Input:    error code and optional int value
// Output:   none
// ---------------------------------------------------------------------------
{
    std::cerr << "CMatrix: Matrix Name: " << m_szName << ".\n";
    switch (nErrorCode)
    {
        case 1:
            std::cerr << "Matrix:: Memory allocation failure.\n";
        break;
        case 2:
            std::cerr << "Matrix::Row-Column reference is out of bounds.\n";
        break;
        case 3:
            std::cerr << "Matrix::Constructor. Invalid number of rows "
                         "or columns.\n";
        break;
        case 4:
            std::cerr << "Matrix::Incompatible matrices.\n";
        break;
    }
    std::cerr << "Unable to populate matrix.\n";
    exit (1);
}

#endif
