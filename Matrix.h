#pragma once

#include <iostream>
#include <random>
#include <typeinfo>

template<typename T = float>
class Matrix {
    static_assert(std::is_arithmetic<T>::value, "T must be a numeric type");

private:
    T** data;
    std::size_t rows;
    std::size_t cols;

    template<typename otherT>
    void check_same_dim(const Matrix<otherT>& other) const{
        if (!(rows == other.getnumrows() && cols == other.getnumcols())) {
            throw std::invalid_argument("matrices must be of the same dimensions");
        }
    }
    void check_inplace(const bool inplace) const{
        if (!inplace){
            throw std::invalid_argument("Inplace if provided, must be true");
        }
    }

public:
    std::string dtype = typeid(T).name();
    /**
     * Test
    */
    Matrix(const std::size_t numRows, const std::size_t numCols){
        rows = numRows;
        cols = numCols;

        // Allocate memory for the matrix
        data = new T*[rows];
        for (int i = 0; i < rows; ++i) {
            data[i] = new T[cols];
        }
    }
    Matrix(const std::size_t dim) : Matrix(dim,dim){}

    ~Matrix() {
        // Deallocate memory for the matrix
        for (std::size_t i = 0; i < rows; ++i) {
            delete[] data[i];
        }
        delete[] data;
    }
    // Getters
    std::size_t getnumrows() const{
        return rows;
    }
    std::size_t getnumcols() const{
        return cols;
    }

    //Conversion constructor (pending)
    template <typename otherT>
    Matrix(const Matrix<otherT>& other) {
        rows = other.getnumrows();
        cols = other.getnumcols();
        // Allocate memory for the matrix
        data = new T*[rows];
        for (std::size_t i = 0; i < rows; ++i){
            data[i] = new T[cols];
        }
        // Copy elements from the other matrix
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                data[i][j] = static_cast<T>(other[i][j]);
            }
        }
    }

    // Static methods
    static Matrix<T> zeros(const std::size_t numRows, const std::size_t numCols){
        Matrix<T> result(numRows, numCols);
        
        // Initialize elements to zero
        for (std::size_t i = 0; i < numRows; ++i) {
            for (std::size_t j = 0; j < numCols; ++j) {
                result[i][j] = T();// Default initialization for the type T
            }
        }
        return result;
    }
    static Matrix<T> zeros(const std::size_t dim){return zeros(dim,dim);}

    static Matrix<T> eye(const std::size_t numRows, const std::size_t numCols){
        Matrix<T> result = Matrix::zeros(numRows, numCols);
        int min = (numRows < numCols) ? numRows : numCols;
        for (std::size_t i = 0; i < min; ++i) {
            result[i][i] = static_cast<T>(1);
        }
        return result;
    }
    static Matrix<T> eye(const std::size_t dim){return eye(dim,dim);}

    static Matrix<T> rand(const std::size_t numRows, const std::size_t numCols) {
        std::random_device dev;
        std::mt19937 gen(dev());
        std::uniform_real_distribution<> dist(0, 1);

        Matrix<T> result(numRows, numCols);

        for (std::size_t i = 0; i < numRows; ++i) {
            for (std::size_t j = 0; j < numCols; ++j) {
                result[i][j] = dist(gen);
            }
        }
        return result;
    }
    static Matrix<T> rand(const std::size_t dim){return rand(dim,dim);}

    // Methods
    Matrix<T> getRow(const std::size_t row) {
        if (row >= rows) {
            throw std::invalid_argument("row exceeds the size of the matrix");
        } else if (row < 0) {
            throw std::invalid_argument("row must be a positive integer");
        }
        Matrix result(1, cols);
        for (std::size_t i = 0; i < cols; ++i) {
            result[0][i] = data[row][i];
        }
        return result;
    }
    Matrix<T> getCol(const std::size_t col) {
        if (col >= cols) {
            throw std::invalid_argument("column exceeds the size of the matrix");
        } else if (col < 0) {
            throw std::invalid_argument("column must be a positive integer");
        }
        Matrix result(rows, 1);
        for (std::size_t i = 0; i < rows; ++i) {
            result[i][0] = data[i][col];
        }
        return result;
    }
    Matrix<T> transpose() const{
        Matrix<T> result(cols, rows);
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = 0; j < cols; ++j){
                result[j][i] = data[i][j];
            }
        }
        return result;
    }
    Matrix<T> transpose(const bool inplace){
        check_inplace(inplace);
        if (rows != cols){
            throw std::logic_error("Matrix must be square to permorm inplace transposition");
        }
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = 0; j < i; ++j){ // Iteration over lower triangle (without diagonal)
                std::swap(data[i][j], data[j][i]);
            }
        }
        return *this;
    }

    Matrix<typename std::common_type<T,float>::type> ref(){
        using resultType = typename std::common_type<T, float>::type;
        Matrix<resultType> result = *this;
        resultType pivot = result[0][0];
        resultType m;
        for (std::size_t i = 1; i < rows; ++i){
            m = result[i][0] / pivot;
            for (std::size_t j = 0; j < cols; ++j){
                result[i][j] -= result[0][j]*m;
            } 
        }
        return result;
    }

    // Operators
    T* operator[](const std::size_t rowIndex) const{
        return data[rowIndex];
    }
    
    friend std::ostream& operator<<(std::ostream& stream, const Matrix<T>& matrix) {
        for (std::size_t i = 0; i < matrix.rows; ++i) {
            for (std::size_t j = 0; j < matrix.cols; ++j) {
                stream << matrix.data[i][j] << " ";
            }
            stream << "\n";
        }
        return stream;
    }
    
    // Conversion operator
    template <typename otherT>
    operator Matrix<otherT>() const {
        return Matrix<otherT>(*this);
    }

    template<typename otherT>
    Matrix<typename std::common_type<T, otherT>::type> operator+(const otherT& num) const {
        static_assert(std::is_arithmetic<otherT>::value, "num must be of numeric type");
        using resultType = typename std::common_type<T, otherT>::type;
        Matrix<resultType> result(rows, cols);

        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] + num;
            }
        }
        return result;
    }
    template<typename otherT>
    Matrix<typename std::common_type<T, otherT>::type> operator+(const Matrix<otherT>& other) const {
        check_same_dim(other);
        using resultType = typename std::common_type<T, otherT>::type;
        Matrix<resultType> result(rows, cols);

        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] + other[i][j];
            }
        }
        return result;
    }

    template<typename otherT>
    Matrix<typename std::common_type<T, otherT>::type> operator-(const otherT& num) const {
        static_assert(std::is_arithmetic<otherT>::value, "num must be of numeric type");
        using resultType = typename std::common_type<T, otherT>::type;
        Matrix<resultType> result(rows, cols);

        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] - num;
            }
        }
        return result;
    }
    template<typename otherT>
    Matrix<typename std::common_type<T, otherT>::type> operator-(const Matrix<otherT>& other) const {
        check_same_dim(other);
        using resultType = typename std::common_type<T, otherT>::type;
        Matrix<resultType> result(rows, cols);

        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] - other[i][j];
            }
        }
        return result;
    }

    template<typename otherT>
    Matrix<typename std::common_type<T, otherT>::type> operator*(const otherT& num) const {
        static_assert(std::is_arithmetic<otherT>::value, "num must be of numeric type");
        using ResultType = typename std::common_type<T, otherT>::type;
        Matrix<ResultType> result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                result[i][j] = num * data[i][j];
            }
        }
        return result;
    }
    template<typename otherT>
    Matrix<typename std::common_type<T, otherT>::type> operator*(const Matrix<otherT>& other) const {
        if (!(cols == other.getnumrows() && rows == other.getnumcols())) {
            throw std::invalid_argument("matrices must have compatible dimensions");
        }
        using ResultType = typename std::common_type<T, otherT>::type;
        Matrix<ResultType> result(rows, other.getnumcols());
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t k = 0; k < cols; ++k) {
                for (std::size_t j = 0; j < other.getnumcols(); ++j) {
                    result[i][j] += data[i][k] * other[k][j];
                }
            }
        }
        return result;
    }

    template <typename otherT>
    Matrix<T> operator^(const otherT& num) const{
        check_same_dim(*this)
        Matrix<T> result = 
    }

    template<typename otherT>
    Matrix<bool> operator<(const otherT& num) const {
        static_assert(std::is_arithmetic<otherT>::value, "num must be of numeric type");
        Matrix<bool> result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] < num;
            }
        }
        return result;
    }
    template<typename otherT>
    Matrix<bool> operator<(const Matrix<otherT>& other) const {
        check_same_dim(other);
        Matrix<bool> result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] < other[i][j];
            }
        }
        return result;
    }

    template<typename otherT>
    Matrix<bool> operator>(const otherT& num) const {
        static_assert(std::is_arithmetic<otherT>::value, "num must be of numeric type");
        Matrix<bool> result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] > num;
            }
        }
        return result;
    }
    template<typename otherT>
    Matrix<bool> operator>(const Matrix<otherT>& other) const {
        check_same_dim(other);
        Matrix<bool> result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                result[i][j] = data[i][j] > other[i][j];
            }
        }
        return result;
    }

    // Other functions

};

