Main function is seq2ge.m

## Units

Units in these scripts follow Pulseq convention, i.e.,  
```
RF: Hz  
Gradients: Hz/m  
Time: sec
```
except when serializing the 'parent block' sequence representation to file
(writeGEseq.m), where hardware units are used (int16).


