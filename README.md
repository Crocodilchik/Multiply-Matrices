# Multiply-Matrices

## Theory
In mathematics, particularly in linear algebra, matrix multiplication is a binary operation that produces a matrix from two matrices. For matrix multiplication, the number of columns in the first matrix must be equal to the number of rows in the second matrix. The resulting matrix, known as the matrix product, has the number of rows of the first and the number of columns of the second matrix. The product of matrices A and B is denoted as AB.

Matrix multiplication was first described by the French mathematician Jacques Philippe Marie Binet in 1812, to represent the composition of linear maps that are represented by matrices. Matrix multiplication is thus a basic tool of linear algebra, and as such has numerous applications in many areas of mathematics, as well as in applied mathematics, statistics, physics, economics, and engineering. Computing matrix products is a central operation in all computational applications of linear algebra.

## Usual multiply
If A is an m × n matrix and B is an n × p matrix:  
<img src="https://github.com/Crocodilchik/Multiply-Matrices/blob/main/%D0%91%D0%B5%D0%B7%D1%8B%D0%BC%D1%8F%D0%BD%D0%BD%D1%8B%D0%B9.png" alt="image" >
the matrix product C = AB (denoted without multiplication signs or dots) is defined to be the m × p matrix
<img src="https://github.com/Crocodilchik/Multiply-Matrices/blob/main/Pic1.png" alt="image2" >
<img src="https://github.com/Crocodilchik/Multiply-Matrices/blob/main/Pic2.png" alt="image3" >

## Vinograd-Coppersmith algorythm
Matrix multiplication according to Winograd (not according to Coppersmith-Winograd!) Consideration of the results of multiplication of two matrices It is obvious that each element in it is a scalar product corresponding to the composition and column of the original matrices. Such a multiplication allows for preprocessing, allowing some of the work to be done ahead of time.

Consider two vectors V = (v1, v2, v3, v4) and W = (w1, w2, w3, w4). Their dot product is: V • W = v1w1 + v2w2 + v3w3 + v4w4.

This equality can be rewritten as: V • W = (v1 + w2)(v2 + w1) + (v3 + w4)(v4 + w3) - v1v2 - v3v4 - w1w2 - w3w4.

Despite the fact that the second number requires more operations than the first: four multiplications - six, instead of three additions - ten, the term on the right side of the reliability predicts a preliminary estimate: its parts can be calculated in advance and memorized each for the conclusion of the first matrix and for each column of the second, which allows performing for each element only the first two multiplications and the subsequent five additions, as well as an additional two additions.

## Strassen algorythm
Strassen’s 7 calls are as follows:

a * (f - h)
(a + b) * h
(c + d) * e
d * (g - e)
(a + d) * (e + h)
(b - d) * (g + h)
(a - c) * (e + f)
Our new matrix C’s new quadrants

matrix C = |p5+p4-p2+p6    p1+p2   |
           |   p3+p4    p1+p5-p3-p7| 
