package eigen_values.methods

import SLAE_methods_functions.MathFunctions
import SLAE_methods_functions.SLAEAccurateMethods
import java.lang.Math.abs
import java.lang.Math.pow

fun powerMethodReversed(matrix: Array<DoubleArray>, accuracy: Double, defShift: Double = -1000.0) {

    var shift = defShift
    var y = Array<DoubleArray>(matrix.size) { DoubleArray(1) {1.0} }

    var normY = MathFunctions.vectorNorm2(y)

    var z = Array<DoubleArray>(matrix.size) { DoubleArray(1) {1.0} }

    var diff: Double
    var normCoef: Double

    for(i in z.indices) {
        z[i][0] = y[i][0]/normY
    }

    do {
        y = SLAEAccurateMethods.GaussMethod(
            MathFunctions.subtractMatrices(
                matrix,
                MathFunctions.multiplyMatrix(MathFunctions.getNeutralMatrix(matrix.size), shift)
            ),
            z
        )

        var mu = 0.0
        var muIter = 0
        for(i in y.indices) {
           if(abs(y[i][0]) > pow(10.0, -8.0)) {
               mu += z[i][0] / y[i][0];
               muIter ++
           }
        }
        mu /= if(muIter > 0) muIter else 1

        shift += mu

        normY = MathFunctions.vectorNorm2(y)

        var nextZ = Array<DoubleArray>(matrix.size) { DoubleArray(1) {1.0} }
        for(i in z.indices) {
            nextZ[i][0] = y[i][0]/normY
        }

        diff = Math.min(
            MathFunctions.vectorNorm2(MathFunctions.subtractMatrices(nextZ, z)),
            MathFunctions.vectorNorm2(MathFunctions.sumMatrices(nextZ, z)),
        )
        normCoef = Math.max(MathFunctions.vectorNorm2(z), MathFunctions.vectorNorm2(nextZ))

        z = nextZ

    } while(diff > accuracy*normCoef)

    print("Eigen value vector: [")
    for(i in z.indices) {
        print("${z[i][0]} ")
    }
    println("]")
    println("Eigen value: $shift")


}