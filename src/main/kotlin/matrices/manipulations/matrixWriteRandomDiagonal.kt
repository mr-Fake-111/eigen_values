package matrices.manipulations

import java.io.File
import java.io.FileWriter
import java.io.IOException

fun matrixWriteRandomDiagonal(dim: Int = 1, file: File) {
    val writer: FileWriter

    try {
        writer = FileWriter(file)
    } catch (e: IOException) {
        throw Exception("Файл не найден или не может быть открыт...")
    }

    if(dim < 1) throw Exception("Размерность матрицы должна быть строго положительной")

    for(i in 0 until dim) {
        for(j in 0 until dim) {
            writer.write( if(i==j) {(Math.random()*20 - 10).toString() + " "} else {"0.0 "})
        }
        writer.write("\n")
    }

    writer.close()
}