/*=============================================================================

    Matt Rasmussen
    Copyright 2007-2011
    Matrix.h
    2007/6/19

=============================================================================*/

#include <stdio.h>

namespace spidir {

template <class ValueType>
ValueType **allocMatrix(int nrows, int ncols)
{
    ValueType **mat = new ValueType* [nrows];
    for (int i=0; i<nrows; i++)
        mat[i] = new ValueType[ncols];
    return mat;
}

template <class ValueType>
void freeMatrix(ValueType **mat, int nrows)
{
    for (int i=0; i<nrows; i++)
        delete [] mat[i];
    delete [] mat;
}



template <class ValueType>
class Matrix
{
public:
    typedef ValueType** ValuePtrType;
    typedef ValueType* ValueRowType;

    Matrix(int nrows, int ncols) :
        m_nrows(nrows),
        m_ncols(ncols)
    {
        m_data = new ValueRowType[nrows];
        for (int i=0; i<nrows; i++)
            m_data[i] = new ValueType[ncols];
    }
    
    ~Matrix()
    {
        for (int i=0; i<m_nrows; i++)
            delete [] m_data[i];
        delete [] m_data;
    }
    
    void setAll(ValueType val)
    {
        for (int i=0; i<m_nrows; i++)
            for (int j=0; j<m_ncols; j++)
                m_data[i][j] = val;
    }
    
    inline int numRows() { return m_nrows; }
    inline int numCols() { return m_ncols; }
    
    inline ValueType* operator[](const int i)
    {
        return m_data[i];
    }
    
    inline ValueType **getMatrix()
    {
        return m_data;
    }
    
    operator ValuePtrType()
    {
        return m_data;
    }
    
private:
    int m_nrows;
    int m_ncols;
    ValueType **m_data;
};



} // namespace spidir
