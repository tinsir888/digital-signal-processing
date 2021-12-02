# Digital-Signal-Processing
Digital signal processing assignments. In the 5th semester in CC NKU, taught by Assoc. Prof. Li Yue.

## Lab 1 Requirements
### Define a signal sequence, require
1. The length of the sequence is limited (the length is self-determined, for example, n=-3:8)
2. The starting position of the sequence is not 0, for example, n=-3 above
3. Able to read and write any position in the sequence, for example x[n], n=-3:8, let x[3]=5; x[4]=x[4]-2, etc.

### Single sequence basic operation
1. Satisfy the pre- and post-zero padding operations
2. Satisfy sequence delay and early operation
3. Satisfy sequence reversal operation
4. Satisfy sequence stretching and compression operations (up-sampling, down-sampling)
5. Satisfy sequence difference and accumulation operations

### Multi-sequence operation
1. Satisfy the addition operation
2. Satisfy the multiplication operation
3. Satisfy the convolution operation
    - Linear convolution
    - Circular convolution
4. Satisfy sequence similarity alignment operation
    - Similarity comparison of sliding window
    - Normalized similarity comparison
