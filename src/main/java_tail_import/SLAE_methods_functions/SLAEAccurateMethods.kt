package SLAE_methods_functions

import SLAE_methods_functions.MathFunctions.getNeutralMatrix
import SLAE_methods_functions.MathFunctions.multiplyMatrices
import SLAE_methods_functions.MathFunctions.multiplyMatrix
import SLAE_methods_functions.MathFunctions.replaceRows
import SLAE_methods_functions.MathFunctions.reversedGaussMove
import SLAE_methods_functions.MathFunctions.straightGaussMove
import SLAE_methods_functions.MathFunctions.subtractMatrices
import SLAE_methods_functions.MathFunctions.sumMatrices
import SLAE_methods_functions.MathFunctions.transpondMatrix
import SLAE_methods_functions.MathFunctions.vectorNorm2
import java.util.*

object SLAEAccurateMethods {
    @Throws(Exception::class)
    fun GaussMethod(
        A: Array<DoubleArray>,
        b: Array<DoubleArray>
    ): Array<DoubleArray> {
        var b = b
        val n = b.size
        var P: Array<DoubleArray> = getNeutralMatrix(n)
        var M = Array(n) { DoubleArray(n) }
        for (i in 0 until n) {
            for (j in 0 until n) {
                M[i][j] = A[i][j]
            }
        } //M == A
        for (i in 0 until n) {
            var max_i = Math.abs(A[i][i])
            var max_i_index = i
            for (j in i until n) {
                if (Math.abs(A[j][i]) > max_i) {
                    max_i_index = j
                    max_i = Math.abs(A[j][i])
                }
            }
            M = replaceRows(M, i, max_i_index)
            P = replaceRows(P, i, max_i_index)
            for (j in i + 1 until n) {
                if(M[i][i] != 0.0) { //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    M[j][i] = M[j][i] / M[i][i]
                }
                for (k in i + 1 until n) {
                    M[j][k] = M[j][k] - M[j][i] * M[i][k]
                }
            }
        }
        val L = getNeutralMatrix(n)
        for (i in 0 until n) {
            for (j in 0 until i) {
                L[i][j] = M[i][j]
            }
        }
        val U =
            subtractMatrices(sumMatrices(M, getNeutralMatrix(n)), L)
        b = multiplyMatrices(P, b)
        val y = reversedGaussMove(L, b)
        return straightGaussMove(U, y)
    }

    @Throws(Exception::class)
    fun ReflectionMethod(
        A: Array<DoubleArray>,
        b: Array<DoubleArray>
    ): Array<DoubleArray> {
        var b = b
        val n = A.size
        var Q: Array<DoubleArray> = getNeutralMatrix(n)
        val R = Array(n) { DoubleArray(n) }

        for (i in 0 until n) {
            R[i] = Arrays.copyOf(A[i], n)
        }

        for (i in 0 until n - 1) {

            val z = Array(n) { DoubleArray(1) }

            for (j in 0 until n) {
                z[j][0] = 0.0
            }
            z[i][0] = 1.0

            val y = Array(n) { DoubleArray(1) }
            for (j in 0 until i) {
                y[j][0] = 0.0
            }

            for (j in i until n) {
                y[j][0] = R[j][i]
            }

            var w = subtractMatrices(y, multiplyMatrix(z, vectorNorm2(y)))
            w = multiplyMatrix(w, 1.0 / vectorNorm2(w))

            val transpondW = Array(1) { DoubleArray(n) }
            for (j in 0 until n) {
                transpondW[0][j] = w[j][0]
            }

            val Q_i = subtractMatrices(
                getNeutralMatrix(n),
                multiplyMatrix(multiplyMatrices(w, transpondW), 2.0)
            )

            var R_i = Array(n) { DoubleArray(n) }
            /*for (j in 0 until n) {
                for (k in 0 until n) {
                    R_i[i][i] = 0.0
                }
            }
            for (j in i until n) {
                for (k in i until n) {
                    R_i[j][k] = R[j][k]
                }
            }*/
            R_i = multiplyMatrices(Q_i, R)

            for (j in i until n) {
                for (k in i until n) {
                    R[j][k] = R_i[j][k]
                }
            }
            Q = multiplyMatrices(Q, Q_i)
        }

        b = multiplyMatrices(transpondMatrix(Q), b)
        return straightGaussMove(R, b)
    }
}