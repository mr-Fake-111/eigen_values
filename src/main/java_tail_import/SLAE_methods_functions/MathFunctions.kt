package SLAE_methods_functions

import java.util.*

object MathFunctions {
    private fun sqrt(`val`: Double, prevRes: Double): Double {
        val newRes = 0.5 * (prevRes + `val` / prevRes)
        return if (Math.abs(newRes - prevRes) < Math.pow(10.0, -9.0)) {
            newRes
        } else {
            sqrt(`val`, newRes)
        }
    }

    fun sqrt(`val`: Double): Double {
        return if (`val` != 1.0) {
            sqrt(`val`, `val`)
        } else 1.0
    }

    fun isDiagonalMajority(matrix: Array<DoubleArray>): Double {
        val n = matrix.size
        var minMajority = Double.MAX_VALUE
        for (i in 0 until n) {
            var `val` = Math.abs(matrix[i][i])
            for (j in 0 until n) {
                if (j != i) {
                    `val` -= Math.abs(matrix[i][j])
                }
            }
            if (`val` <= 0) return (-1).toDouble()
            if (`val` < minMajority) minMajority = `val`
        }
        return minMajority
    }

    @JvmStatic
    fun transpondMatrix(matrix: Array<DoubleArray>): Array<DoubleArray> {
        val n = matrix.size
        val transpondedMatrix = Array(n) { DoubleArray(n) }
        for (i in 0 until n) {
            for (j in 0 until n) {
                transpondedMatrix[i][j] = matrix[j][i]
            }
        }
        return transpondedMatrix
    }

    fun norm1(matrix: Array<DoubleArray>): Double {
        var maxColumnValue = 0.0
        for (j in matrix[0].indices) {
            var columnSum = 0.0
            for (i in matrix.indices) {
                columnSum += Math.abs(matrix[i][j])
            }
            if (columnSum > maxColumnValue) maxColumnValue = columnSum
        }
        return maxColumnValue
    }

    fun normInfinity(matrix: Array<DoubleArray>): Double {
        val n = matrix.size
        return if (matrix[0].size > 1) {
            var maxRowValue = 0.0
            for (i in 0 until n) {
                var rowSum = 0.0
                for (j in 0 until n) {
                    rowSum += Math.abs(matrix[i][j])
                }
                if (rowSum > maxRowValue) maxRowValue = rowSum
            }
            maxRowValue
        } else {
            var maxValue = 0.0
            for (i in 0 until n) {
                if (maxValue < Math.abs(matrix[i][0])) maxValue = Math.abs(matrix[i][0])
            }
            maxValue
        }
    }

    @JvmStatic
    @Throws(Exception::class)
    fun sumMatrices(m1: Array<DoubleArray>, m2: Array<DoubleArray>): Array<DoubleArray> {
        if (m1.size != m2.size || m1[0].size != m2[0].size) throw Exception("matrices are not compatible for sum")
        val resMatrix = Array(m1.size) { DoubleArray(m1[0].size) }
        for (i in m1.indices) {
            for (j in m1[0].indices) {
                resMatrix[i][j] = m1[i][j] + m2[i][j]
            }
        }
        return resMatrix
    }

    @JvmStatic
    @Throws(Exception::class)
    fun subtractMatrices(m1: Array<DoubleArray>, m2: Array<DoubleArray>): Array<DoubleArray> {
        return sumMatrices(m1, multiplyMatrix(m2, -1.0))
    }

    @JvmStatic
    @Throws(Exception::class)
    fun multiplyMatrices(m1: Array<DoubleArray>, m2: Array<DoubleArray>): Array<DoubleArray> {
        if (m1[0].size != m2.size) throw Exception("matrices are not compatible for multiplication")
        val resMatrix = Array(m1.size) { DoubleArray(m2[0].size) }
        for (i in resMatrix.indices) {
            for (j in resMatrix[0].indices) {
                var value = 0.0
                for (k in m2.indices) {
                    value += m1[i][k] * m2[k][j]
                }
                resMatrix[i][j] = value
            }
        }
        return resMatrix
    }

    @JvmStatic
    fun multiplyMatrix(m: Array<DoubleArray>, value: Double): Array<DoubleArray> {
        val resMatrix = Array(m.size) { DoubleArray(m[0].size) }
        for (i in m.indices) {
            for (j in m[0].indices) {
                resMatrix[i][j] = value * m[i][j]
            }
        }
        return resMatrix
    }

    @JvmStatic
    fun getNeutralMatrix(n: Int): Array<DoubleArray> {
        val matrix = Array(n) { DoubleArray(n) }
        for (i in 0 until n) {
            for (j in 0 until n) {
                if (i == j) matrix[i][j] = 1.0 else matrix[i][j] = 0.0
            }
        }
        return matrix
    }

    fun getBadMatrix(n: Int, epsilon: Double): Array<DoubleArray> {
        val badMatrix = getNeutralMatrix(n)
        for (i in 0 until n) {
            for (j in i + 1 until n) {
                badMatrix[i][j] = -1 - epsilon * 10
            }
        }
        for (i in 0 until n) {
            badMatrix[i][i] += epsilon * 10
        }
        for (i in 0 until n) {
            for (j in 0 until i) {
                badMatrix[i][j] = epsilon * 10
            }
        }
        return badMatrix
    }

    fun printMatrix(matrix: Array<DoubleArray>) {
        for (i in matrix.indices) {
            for (j in matrix[0].indices) {
                print(matrix[i][j].toString() + " ")
            }
            println()
        }
    }

    @JvmStatic
    fun straightGaussMove(matrix: Array<DoubleArray>, b: Array<DoubleArray>): Array<DoubleArray> {
        val x = Array(b.size) { DoubleArray(1) }
        val n = x.size
        for (i in n - 1 downTo 0) {
            x[i][0] = b[i][0]
            for (j in i + 1 until n) {
                x[i][0] -= x[j][0] * matrix[i][j]
            }
            if(matrix[i][i] != 0.0) {
                x[i][0] /= matrix[i][i]
            }
        }
        return x
    }

    @JvmStatic
    fun reversedGaussMove(matrix: Array<DoubleArray>, b: Array<DoubleArray>): Array<DoubleArray> {
        val x = Array(b.size) { DoubleArray(1) }
        val n = x.size
        for (i in 0 until n) {
            x[i][0] = b[i][0]
            for (j in i - 1 downTo 0) {
                x[i][0] -= x[j][0] * matrix[i][j]
            }
            if(matrix[i][i] != 0.0) {
                x[i][0] /= matrix[i][i]
            }
        }
        return x
    }

    @JvmStatic
    fun replaceRows(matrix: Array<DoubleArray>, i: Int, j: Int): Array<DoubleArray> {
        val help = Arrays.copyOf(matrix[i], matrix[i].size)
        matrix[i] = matrix[j]
        matrix[j] = help
        return matrix
    }

    @JvmStatic
    @Throws(Exception::class)
    fun vectorNorm2(vector: Array<DoubleArray>): Double {
        if (vector[0].size != 1) throw Exception("Matrix is not a vector-column")
        var res = 0.0
        for (i in vector.indices) {
            res += Math.pow(vector[i][0], 2.0)
        }
        res = Math.sqrt(res)
        return res
    }
}