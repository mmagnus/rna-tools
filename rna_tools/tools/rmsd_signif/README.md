rnastruc_p_value
-------------------------------------------------------------------------------

A wrapper for RNA_PredictionSignificance by:

> Hajdin, C. E., Ding, F., Dokholyan, N. V, & Weeks, K. M. (2010). On the significance of an RNA tertiary structure prediction. RNA (New York, N.Y.), 16(7), 1340–9. doi:10.1261/rna.1837410

Test the code:

    ➜  rnastruc_p_value git:(master) ✗ ./rnastruc_p_value.py
    Without base pair constraints
     The average RMSD by chance: 97.6 
     Z-score of the prediction: -53.67
     p-value of the prediction: < 1.00e-06
    
    With base pair constraints
     The average RMSD by chance: 70.6 
     Z-score of the prediction: -38.68
     p-value of the prediction: < 1.00e-06
    
    [1e-06, 1e-06]

Compile the source code

    g++ RNA_PredictionSignificance.cpp -o RNA_PredictionSignificance.app # for Mac & Linux



