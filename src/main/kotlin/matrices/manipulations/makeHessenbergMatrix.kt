package matrices.manipulations

import SLAE_methods_functions.MathFunctions
import java.lang.Math.pow
import java.util.Arrays
import java.util.DoubleSummaryStatistics
import kotlin.math.sign
import kotlin.math.sqrt

fun makeHessenbergMatrix(A: Array<DoubleArray>): Array<DoubleArray> {

    val n = A.size
    var B = Array<DoubleArray>(n) {DoubleArray(n)}
    for (i in B.indices) {
        B[i] = Arrays.copyOf(A[i], n)
    }

    for( i in 1 until n-1) {
        var v = Array<DoubleArray>(n) { doubleArrayOf(0.0) }
        var s = 0.0
        for(j in i until n) {
            v[j][0] = B[j][i-1]
            s += pow(B[j][i-1], 2.0)
        }
        s = sign(B[i][i-1])*sqrt(s)

        v[i][0] -= s
        v = MathFunctions.multiplyMatrix(v, 1/sqrt(2*s*(s-B[i][i-1])))

        val transpondV = Array(1) { DoubleArray(n) }
        for (j in 0 until n) {
            transpondV[0][j] = v[j][0]
        }

        val H = MathFunctions.subtractMatrices(
            MathFunctions.getNeutralMatrix(n),
            MathFunctions.multiplyMatrix(MathFunctions.multiplyMatrices(v, transpondV), 2.0)
        )

        B = MathFunctions.multiplyMatrices(H, B)
        B = MathFunctions.multiplyMatrices(B, H)
    }

    return B
}