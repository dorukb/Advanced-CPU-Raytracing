#ifndef __DORKTRACER_MATRIX__
#define __DORKTRACER_MATRIX__

#include <vector>
#include <stdint.h>
#include <iostream>
#include "helperMath.h"

namespace DorkTracer{

    typedef std::vector<std::vector<double>> mat;

    class Matrix{
        public:
            int row, col;

            Matrix(){
                row = col = 0;
            }

            Matrix(int row, int col)
            {
                m_mat.resize(row, std::vector<double>(col, 0.0f));
                this->row = row;
                this->col = col;
            } 
            
            static Matrix GetTranslation(double tx, double ty, double tz){
                Matrix t(4,4);
                t[0][0] = t[1][1] = t[2][2] = t[3][3] = 1.0f;
                t[0][3] = tx;
                t[1][3] = ty;
                t[2][3] = tz;
                return t;
            }

            static Matrix GetScale(double sx, double sy, double sz){
                Matrix s(4,4);
                s[0][0] = sx;
                s[1][1] = sy;
                s[2][2] = sz;
                s[3][3] = 1.0f;
                return s;
            }

            static Matrix GetRotationAroundX(double angle)
            {
                Matrix r(4,4);
                r[0][0] = r[3][3] = 1.0f;
                r[1][1] = r[2][2] = std::cos(angle);
                r[1][2] = -std::sin(angle);
                r[2][1] = std::sin(angle);
                return r;
            }

            static Matrix GetRotationAroundY(double angle)
            {
                Matrix rx(4,4);
                rx[0][0] = rx[2][2] = std::cos(angle);
                rx[1][1] = rx[3][3] = 1.0f;
                rx[0][2] = std::sin(angle);
                rx[2][0] = -std::sin(angle);
                return rx;
            }

            static Matrix GetRotationAroundZ(double angle)
            {
                Matrix r(4,4);
                r[0][0] = r[1][1] = std::cos(angle);
                r[2][2] = r[3][3] = 1.0f;
                r[1][0] = std::sin(angle);
                r[0][1] = -std::sin(angle);
                return r;
            }
            
            static Matrix GetTranspose(Matrix& m){
                Matrix result(m.col, m.row);
                for(int i = 0; i < m.row; i++){
                    for(int j = 0; j < m.col; j++){
                        result[j][i] = m[i][j];
                    }
                }
                return result;
            }   

            std::vector<double>& operator[](uint32_t idx) {
                return m_mat.at(idx);
            }

            static Vec3f ApplyTransform(Matrix& t, Vec4f v){
                // only supports 4x4 matrices.
                if(t.col != 4 && t.row != 4){
                    std::cout << "Cannot apply transform. given t is not 4x4." << std::endl;
                    return Vec3f();
                }

                Vec3f res;
                res.x = t[0][0] * v.x + t[0][1] * v.y +
                        t[0][2] * v.z + t[0][3] * v.w;

                res.y = t[1][0] * v.x + t[1][1] * v.y +
                        t[1][2] * v.z + t[1][3] * v.w;

                res.z = t[2][0] * v.x + t[2][1] * v.y +
                        t[2][2] * v.z + t[2][3] * v.w;

                float w = t[3][0] * v.x + t[3][1] * v.y +
                        t[3][2] * v.z + t[3][3] * v.w;

                return res;
            }
            
            static Vec3f ApplyTransformToPoint(Matrix& t, Vec3f p)
            {
                return ApplyTransform(t, Vec4f(p, 1.0f));
            }

            static Vec3f ApplyTransformToVector(Matrix& t, Vec3f p)
            {
                return ApplyTransform(t, Vec4f(p, 0.0f));
            }

            Matrix operator*(Matrix& rhs){
                Matrix result(this->row, rhs.col);
                if(this->col != rhs.row){
                    std::cout << "Invalid multiplication request, dimensions do not match." << std::endl;
                    return *this;
                }

                for(int i = 0; i < this->row; i++)
                {
                    for(int j = 0 ; j < rhs.col; j++)
                    {
                        result[i][j] = 0.0f;
                        for(int k = 0; k < this->col; k++)
                        {
                            result[i][j] += m_mat[i][k] * rhs[k][j];
                        }
                    }
                }
                return result;
            }
            void MakeIdentity(){
                // Only works for square matrices.
                if(this->row != this->col) return;

                for(int i = 0; i < this->row; i++)
                {
                    m_mat[i][i] = 1.0f;
                }
            }

        private:
            mat m_mat;
    };
}

#endif