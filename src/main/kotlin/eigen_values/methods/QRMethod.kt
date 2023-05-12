package eigen_values.methods

import SLAE_methods_functions.MathFunctions
import matrices.manipulations.makeHessenbergMatrix
import java.lang.Math.abs
import java.lang.Math.pow
import java.util.*

fun QRMethod(
    A: Array<DoubleArray>,
    eigens: ArrayList<Double>
): ArrayList<Double> {

    val n = A.size
    var Q: Array<DoubleArray> = MathFunctions.getNeutralMatrix(n)
    val R = Array(n) { DoubleArray(n) }
    var resMat = Array(n) { DoubleArray(n) }
    var iter = 0

    for(i in 0 until n) {
        resMat[i] = Arrays.copyOf(A[i], n)
    }
    var swish: Array<DoubleArray>

    do {
        swish = MathFunctions.multiplyMatrix(MathFunctions.getNeutralMatrix(n), resMat[n-1][n-1])

        resMat = MathFunctions.subtractMatrices(resMat, swish)

        for (i in 0 until n) {
            R[i] = Arrays.copyOf(resMat[i], n)
        }
        var Q: Array<DoubleArray> = MathFunctions.getNeutralMatrix(n)

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

            var w = MathFunctions.subtractMatrices(y, MathFunctions.multiplyMatrix(z, MathFunctions.vectorNorm2(y)))
            w = MathFunctions.multiplyMatrix(w, 1.0 / MathFunctions.vectorNorm2(w))

            val transpondW = Array(1) { DoubleArray(n) }
            for (j in 0 until n) {
                transpondW[0][j] = w[j][0]
            }

            val Q_i = MathFunctions.subtractMatrices(
                MathFunctions.getNeutralMatrix(n),
                MathFunctions.multiplyMatrix(MathFunctions.multiplyMatrices(w, transpondW), 2.0)
            )

            var R_i = Array(n) { DoubleArray(n) }
            R_i = MathFunctions.multiplyMatrices(Q_i, R)

            for (j in i until n) {
                for (k in i until n) {
                    R[j][k] = R_i[j][k]
                }
            }
            Q = MathFunctions.multiplyMatrices(Q, Q_i)
        }

        resMat = MathFunctions.multiplyMatrices(R, Q)
        resMat = MathFunctions.sumMatrices(resMat, swish)
        iter ++

    } while(abs(resMat[n-1][n-2]) > pow(10.0, -8.0) && iter < 1000)

    eigens.add(resMat[n-1][n-1])

    if(n > 2) {
        val lowA = Array<DoubleArray>(n - 1) { DoubleArray(n - 1) }
        for (i in 0 until n-1) {
            for (j in 0 until n-1) {
                lowA[i][j] = resMat[i][j]
            }
        }
        QRMethod(lowA, eigens)
    } else {
        eigens.add(resMat[0][0])
    }

    return eigens
}