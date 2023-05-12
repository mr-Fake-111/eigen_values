import SLAE_methods_functions.SLAEAccurateMethods
import eigen_values.methods.QRMethod
import eigen_values.methods.powerMethod
import eigen_values.methods.powerMethodReversed
import matrices.manipulations.*
import java.io.File
import java.lang.Math.pow
import java.util.*

fun main(args:Array<String>) {

    val scan = Scanner(File("src/main/resources/matrixC.txt"))

    var matrix = Array<DoubleArray>(7) {DoubleArray(7) {0.0} }
    for (i in 0..6) {
        for (j in 0..6) {
            matrix[i][j] = scan.next().toDouble()
        }
    }
    matrix = makeHessenbergMatrix(matrix)

    powerMethod(matrix, pow(10.0, -10.0))
    powerMethodReversed(matrix, pow(10.0, -10.0), 23.0)


    val eigens = ArrayList<Double>()
    QRMethod(matrix, eigens)
    for(i in eigens) {
        print(i)
        print(" ")
    }
}