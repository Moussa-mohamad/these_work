
import os
import pandas as pd  # for pandas.DataFrame objects (https://pandas.pydata.org/)
import numpy as np  # for numpy.matrix objects (https://numpy.org/)
from amplpy import AMPL
from scipy.sparse import csr_matrix

def kin_ipopt( data,A,b,c,length_edge):
    # You can install amplpy with "python -m pip install amplpy"
    

    # Create an AMPL instance
    ampl = AMPL()
   
    # Set the solver to use
    solver = "highs" 
    ampl.set_option("solver", "ipopt")
    
    ampl.eval(
        r"""
        set Aline;
        set Gline;
        set bline;
        set hline;
        set cline;
        set lenline;


        param n := 3 ;

        param startTime >=0;
        param endTime >=0;
        param sumcons ;
        param ind ;
        param initind;
        param k;
        param count;

        param AColInd {Aline} >= 0;
        param ALineInd {Aline} >= 0;
        param AVal {Aline} ;

        param LenVal {lenline} ;


        param BVal {bline} ;


        param A{Aline,1..n} ;


        

        param len{lenline,1..5} ;

        param numbZ ;
        param FirstIndd  ;
        param FirstIndZ  ;
        param m :=5;
        param avail {Aline} ;



        


        param c{cline,1..2} ;

        param CInd{cline} >=0;
        param CVal{cline} ;



        
    """
    )
    
 
    

    # 1. Send the data from "A" to AMPL and initialize the indexing set "Aline"
    ampl.set_data(A, "Aline")
    # 2. Send the data from "b" to AMPL and initialize the indexing set "bline"
    ampl.set_data(b, "bline")
    # 3. Send the data from "c" to AMPL and initialize the indexing set "cline"
    ampl.set_data(c, "cline")
    # 4. Send the data from "length_edge" to AMPL and initialize the indexing set "lenline"
    ampl.set_data(length_edge, "lenline")
    
    ampl.eval(
        r"""
        
       param A_lnum:= max {i in Aline} ALineInd[i];

        param A_cnum:= max {i in Aline} AColInd[i];

       
        var gamma1{1.. numbZ/6} >=0;
        var gamma2{1.. numbZ/6} >=0;


        let sumcons := 0;
        let ind :=1;
        
        param nb = (A_lnum-1)/3 - numbZ/2;
        
        let numbZ := A_cnum-A_lnum+1;


       var x_cin {1..FirstIndZ + numbZ -1};
        var x{1..nb*3 + numbZ};
        var py{1..numbZ/6};

        let FirstIndd := A_lnum+1 ;
        let FirstIndZ :=(A_cnum+1-A_lnum)/2+2 ;



        param nc = 1000;

        param fc = 500;
        param fcsol = 95;
        param cp = 1;
        param ft = 0;
        param e = 0;

        param allinitind{1..A_lnum};
        param allind{1..A_lnum};



        for {i in 1..A_lnum} {
        let initind := ind; 
        repeat while ALineInd[ind] == i   { if ind+1 > card(Aline) then break; else {let ind := ind +1;} 
          } 
          if ind+1 > card(Aline) then let allind[i] := ind; else let allind[i] := ind-1;
          let allinitind[i] := initind;

          }


subject to equalities:
sum {j in allinitind[A_lnum-numbZ/2]..allind[A_lnum-numbZ/2]} AVal[j]*x[AColInd[j]+1  -   (FirstIndZ + numbZ-1) ]   = BVal[A_lnum-numbZ/2];


subject to equalities1 {i in A_lnum-numbZ/2+1..A_lnum} :
if ((i - A_lnum+numbZ/2-1) mod 3) ==0  then sum {j in allinitind[i]..allind[i]} AVal[j]*x[AColInd[j]+1 - (FirstIndZ + numbZ-1)  ]   = BVal[i];

subject to equalities2 {i in A_lnum-numbZ/2+1..A_lnum} :
if ((i - A_lnum+numbZ/2-1) mod 3) ==1  then sum {j in allinitind[i]..allind[i]} AVal[j]*x[AColInd[j]+1 -  (FirstIndZ + numbZ-1)] + (gamma1[(i - A_lnum+numbZ/2+1)/3] + gamma2[(i - A_lnum+numbZ/2+1)/3])*(-(fc-ft)*LenVal[(i - A_lnum+numbZ/2+1)/3]/(2*(ft+fc)) + py[(i - A_lnum + numbZ/2-1 + 2)/3]/(fc+ft))    = BVal[i];

subject to equalities3 {i in A_lnum-numbZ/2+1..A_lnum} :
if ((i - A_lnum+numbZ/2-1) mod 3) ==2  then sum {j in allinitind[i]..allind[i]} AVal[j]*x[AColInd[j]+1 -  (FirstIndZ + numbZ-1) ] + (gamma1[(i - A_lnum+numbZ/2)/3] - gamma2[(i - A_lnum+numbZ/2)/3])*LenVal[(i - A_lnum+numbZ/2)/3]/2   = BVal[i];



subject to inequalities1 {i in 0..numbZ-1}:
x[FirstIndd+i - (FirstIndZ + numbZ-1)] >= 0;


minimize total_cost: 
sum {i in 1..3*nb} x[i+FirstIndZ + numbZ-1 - (FirstIndZ + numbZ-1) ]*BVal[i] + sum {i in 1..numbZ} x[FirstIndZ + numbZ + 3*nb -1 + i  - (FirstIndZ + numbZ-1)]*BVal[3*nb + i] + sum {i in 1.. numbZ/6} py[i]*py[i]*(gamma1[i]+ gamma2[i])/(2*(ft+fc))  + sum {i in 1.. numbZ/6} (gamma1[i]+ gamma2[i])*0.5*ft*fc*LenVal[i]*LenVal[i]/(ft+fc) ;

    """
    )

   
    
    # Solve
    ampl.solve()

    
    # Get objective entity by AMPL name
    obj = ampl.get_objective("total_cost").value()
    result = ampl.getValue("solve_result") 
        
   
    # Print the objective
    print("Objective is:", obj)

    x = ampl.get_variable("x")
    nb = int(ampl.get_parameter("nb").value() )
    numbZ = int(ampl.get_parameter("numbZ").value())
    nc = int(numbZ/6)
    
    statsol = np.empty((0,nc*3+1));
 
    # Print it
    for i in range(0,3*nc+1):
        statsol = np.append(statsol,x[i+1].value()) 

    
    return statsol,x,nb,nc,obj,result
