package eigen_values.methods

import SLAE_methods_functions.MathFunctions
import java.lang.Math.*
import java.util.*
import java.util.Arrays.copyOf

//не работает для диагональных и сильно разряженных матриц (так как данный случай приводит к неизменяемости вектора y,
//и, как следствие, стагнации алгоритма). Т. е. алгоритм завершается, однако координаты вектора lambda не задают усреднение одного максимального числа,
//а содержат все собственные числа разом.
//в связи с этим вычисление собственного вектора также невозможно
//
//Для нетривиальных матриц все работает корректно
fun powerMethod(matrix: Array<DoubleArray>, accuracy: Double): Double {

    var globalIter = 0

    var y = Array<DoubleArray>(matrix.size) { DoubleArray(1) {1.0} }

    var normY = MathFunctions.vectorNorm2(y)

    var z = Array<DoubleArray>(matrix.size) { DoubleArray(1) {1.0} }
    var lambda = Array<DoubleArray>(matrix.size) { DoubleArray(1) {0.0} }

    for(i in z.indices) {
        z[i][0] = y[i][0]/normY
    }
    y = MathFunctions.multiplyMatrices(matrix, z)

    for(i in lambda.indices) {
        if(abs(z[i][0]) > pow(10.0, -8.0)) {
            lambda[i][0] = y[i][0] / z[i][0]
        }
    }

    var newLambda = Array<DoubleArray>(matrix.size) { DoubleArray(1) {0.0} }

    var diff: Double
    var normCoef: Double

    do {
        normY = MathFunctions.vectorNorm2(y)

        for(i in z.indices) {
            z[i][0] = y[i][0]/normY
        }
        y = MathFunctions.multiplyMatrices(matrix, z)


        for(i in newLambda.indices) {
            if(abs(z[i][0]) > pow(10.0, -8.0)) {
                newLambda[i][0] = y[i][0] / z[i][0]
            }
        }

        diff = MathFunctions.vectorNorm2(MathFunctions.subtractMatrices(newLambda, lambda))
        normCoef = max(MathFunctions.vectorNorm2(lambda), MathFunctions.vectorNorm2(newLambda))

        for(i in lambda.indices) {
            lambda[i][0] = newLambda[i][0]
        }

        globalIter ++

    } while(diff > accuracy*normCoef)

    var maxEigenValue = Double.MIN_VALUE
    var iter  = 0
    for(i in newLambda.indices) {
        if(newLambda[i][0] > maxEigenValue) {
            maxEigenValue = newLambda[i][0]
        }
    }

    val eigenVector = Array<DoubleArray>(matrix.size) { DoubleArray(1) {0.0} }
    for(i in eigenVector.indices) {
        eigenVector[i][0] = y[i][0]/pow(maxEigenValue, globalIter.toDouble())
    }
    val normEigenVector = MathFunctions.vectorNorm2(eigenVector)

    print("max eigen value vector: [")
    for(i in eigenVector.indices) {
        print("${eigenVector[i][0]/normEigenVector} ")
    }
    println("]")
    println("Max eigen value: $maxEigenValue")

    return maxEigenValue
}