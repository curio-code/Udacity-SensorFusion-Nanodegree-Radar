# SFND Radar Final Project

## Implemenation of 2D-CFAR
1. Looping through every Cell under test (CUT)
```
for i   = T_h+G_h+1:Nr/2-T_h-G_h
    for j = T_v+G_v+1:Nd-T_v-G_v
```
  * here ```T_h```, ```T_v``` are the training cells ```G_h```, ```G_v``` are gaurds cells around the CUT.
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
