# Sequence Encodings for Predictive Monitoring
Supplementary material for the article "Complex Symbolic Sequence Encodings for Predictive Monitoring of Business Processes" by Anna Leontjeva, Raffaele Conforti, Chiara Di Francescomarino, Marlon Dumas, Fabrizio Maria Maggi (http://link.springer.com/chapter/10.1007/978-3-319-23063-4_21).

This paper addresses the problem of predicting the outcome of an ongoing case of a business process based on event logs. We explore different encoding schemes for the complex sequences so that it can be translated into the feature space suitable for the supervised models.  Different combinations of encodings have been tried and compared: 

* boolean
* frequency-based 
* simple index encoding
* index latest payload encoding
* HMM-based encoding.

Random forest is trained on the features extracted from all of the encodings. 

We also provide one example dataset based on a BPIC 2011 data (http://www.win.tue.nl/bpi/2011/challenge). Note that this code is a supplementary material and should be modified for the use with other datasets. 
If you use the code please cite the following paper:
```
@incollection{leontjeva2015complex,
  title={Complex Symbolic Sequence Encodings for Predictive Monitoring of Business Processes},
  author={Leontjeva, Anna and Conforti, Raffaele and Di Francescomarino, Chiara and Dumas, Marlon and Maggi, Fabrizio Maria},
  booktitle={Business Process Management},
  pages={297--313},
  year={2015},
  publisher={Springer}
}
```
