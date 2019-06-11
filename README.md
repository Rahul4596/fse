# fse
My reimplementation of Yann Collet's FSE 

This is a FSE coder based on Yann Collet's blog. Currently, only grayscale images can be compressed. This implementation is slower than the reference and as efficient. FOr higher file sizes, this implementation is slower but more efficient. NOw, the whole stream is divided into blocks.

Things to try:
completely avoiding multiplication or division
minimising loops in beat read and write
track 2 states simultaneously

