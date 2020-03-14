# SFND Radar Final Project
## Project Outline
![alt text](https://github.com/curio-code/Udacity-SensorFusion-Nanodegree-Radar/blob/master/media/projectOutline.png)

Following task were completed in this project:-
1. Configured the FMCW waveform based on the system requirements.
2. Defined the range and velocity of target and simulate its displacement.
3. For the same simulation loop processed the transmit and receive signal to determine the beat signal
4. Performed Range FFT on the received signal to determine the Range
5. Towards the end, performed the CFAR processing on the output of 2nd FFT to display the target.

## 2D - FFT
MATLAB inbuilt funtion was utlilized for this purpose ```sig_fft2 = fft2(Mix,Nr,Nd)```
![alt text](https://github.com/curio-code/Udacity-SensorFusion-Nanodegree-Radar/blob/master/media/2dfft.png)

## Implemenation of 2D-CFAR
1. Looping through every Cell under test (CUT)
```
for i   = T_h+G_h+1:Nr/2-T_h-G_h
    for j = T_v+G_v+1:Nd-T_v-G_v
```
  * here ```T_h```, ```T_v``` are the training cells ```G_h```, ```G_v``` are gaurds cells around the CUT.
![alt text](https://github.com/curio-code/Udacity-SensorFusion-Nanodegree-Radar/blob/master/media/2dcfar_cell.png)

2. Calculating the noise in the training cells around the CUT.
```
S1 = sum(db2pow(RDM(i-T_h-G_h:i+T_h+G_h, j-T_v-G_v:j+T_v+G_v)),[1 2]);    
S2 = sum(db2pow(RDM(i-G_h:i+G_h, j-G_v:j+G_v)),[1 2]);
train_total = S1-S2;
```    
3. Now, averaging the total noise ``` train_total``` in previous step and adding a offset to that determine and a threshold and supress the noises at the edges.
 ```
 threshold = train_total/train_cells;
threshold = pow2db(threshold) + offset;
threshold = db2pow(threshold);
```
4. And finally comparing with the threshold with the CUT.
```
if (signal <= threshold)
            signal = 0;
        else 
            signal = 1;
end
```
## Selection of Training, Guard cells and Offset
  * Training cells selected across the rows and columns: ```T_h```, ```T_v``` are **12** and **6**.
  * Gaurd cells selected across the rows and columns: ```G_h```, ```G_v``` are **6** and **3**.
  * Offset was set to 5 but selecting was an iterative process, it was compared with previous results for most optimal value.
    * offset = **5**
   
## Suppression of the non-thresholded cells at the edges
For the final output of CFAR ```signal_cfar``` I inherently considered a zeros matrix of the dimension of RDM ```zeros(Nr/2,Nd)```, so in the end I don't have to consider the cells at the edges.
