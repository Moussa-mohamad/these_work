import os
import pandas as pd  # for pandas.DataFrame objects (https://pandas.pydata.org/)
import numpy as np  # for numpy.matrix objects (https://numpy.org/)
from amplpy import AMPL

def stat_ampl( A,b,c,length_edge):
    # You can install amplpy with "python -m pip install amplpy"
    

    # Create an AMPL instance
    ampl = AMPL()
   
    # Set the solver to use
    solver = "highs" 
    ampl.set_option("solver", "ipopt")
    
    ampl.eval(
        r"""
        set Aline;
        set bline;
        set cline;
        set lenline;


        param n := 3 ;

        param startTime >=0;
        param endTime >=0;
        param sumcons ;
        param ind ;
        param initind;
        param k;


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



        var x{1..FirstIndZ + numbZ -1};



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


        let ind :=1;
        
        param nb = (A_lnum-1)/3 - numbZ/2;
        
        let numbZ := A_cnum-A_lnum+1;

        let FirstIndd := A_lnum+1 ;
        let FirstIndZ :=(A_cnum+1-A_lnum)/2+2 ;

        var alpha1{1..numbZ/6};
        var alpha2{1..numbZ/6};



        param nc = 1000;

        param fc = 500;
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



        subject to equalities1 {i in 1..nb*3} :
        sum {j in allinitind[i]..allind[i]} AVal[j]*x[AColInd[j]+1] = BVal[i]; 

        #subject to equalities2 {i in nb*3+1..A_lnum-numbZ/2-1} :
        #if (((i-nb*3-1) mod 6) <> 4 and  ((i-nb*3-1) mod 6) <> 5   and (i-nb*3-1 < 6*0 or i-nb*3-1 >6*0+5 ) and (i-nb*3-1 < 6*4 or i-nb*3-1 >6*4+5 ) and (i-nb*3-1 < 6*8 or i-nb*3-1 >6*8+5 ) and (i-nb*3-1 < 6*11 or i-nb*3-1 >6*11+5 ) and (i-nb*3-1 < 6*15 or i-nb*3-1 >6*15+5 ) and (i-nb*3-1 < 6*19 or i-nb*3-1 >6*19+5 ) and (i-nb*3-1 < 6*23 or i-nb*3-1 >6*23+5 ) and (i-nb*3-1 < 6*25 or i-nb*3-1 >6*25+5 ) and (i-nb*3-1 < 6*27 or i-nb*3-1 >6*27+5 ) and (i-nb*3-1 < 6*29 or i-nb*3-1 >6*29+5 )) then 
        #sum {j in allinitind[i]..allind[i]} AVal[j]*x[AColInd[j]+1] = if (((i-nb*3-1) mod 6) <> 4 and  ((i-nb*3-1) mod 6) <> 5  and (i-nb*3-1 < 6*0 or i-nb*3-1 >6*0+5 ) and (i-nb*3-1 < 6*4 or i-nb*3-1 >6*4+5 ) and (i-nb*3-1 < 6*8 or i-nb*3-1 >6*8+5 ) and (i-nb*3-1 < 6*11 or i-nb*3-1 >6*11+5 ) and (i-nb*3-1 < 6*15 or i-nb*3-1 >6*15+5 ) and (i-nb*3-1 < 6*19 or i-nb*3-1 >6*19+5 ) and (i-nb*3-1 < 6*23 or i-nb*3-1 >6*23+5 ) and (i-nb*3-1 < 6*25 or i-nb*3-1 >6*25+5 ) and (i-nb*3-1 < 6*27 or i-nb*3-1 >6*27+5 ) and (i-nb*3-1 < 6*29 or i-nb*3-1 >6*29+5 )) then BVal[i] else  0; 


        subject to equalities2 {i in nb*3+1..A_lnum-numbZ/2-1} :
        if ( ((i-nb*3-1) mod 6) <> 2 and ((i-nb*3-1) mod 6) <> 3 and ((i-nb*3-1) mod 6) <> 4 and  ((i-nb*3-1) mod 6) <> 5  ) then sum {j in allinitind[i]..allind[i]} AVal[j]*x[AColInd[j]+1] = if ( ((i-nb*3-1) mod 6) <> 2 and ((i-nb*3-1) mod 6) <> 3 and ((i-nb*3-1) mod 6) <> 4 and  ((i-nb*3-1) mod 6) <> 5  ) then BVal[i]  else 0; 

        
        subject to inequalities2 {i in 0..numbZ -1}:
        x[FirstIndZ+i] >= 0;

        
        subject to inequalities3 {i in 1..numbZ/6}:
        alpha1[i] >= 0;

        
        subject to inequalities4 {i in 1..numbZ/6}:
        alpha2[i] >= 0;


        
         
        subject to crushing_top {i in lenline}:
        x[3*i]*LenVal[i]/2 - x[3*i-1]*((fc-ft)*LenVal[i]/(2*(ft+fc)) -x[3*i-1]/(2*(fc+ft) )) - 0.5*fc*ft*LenVal[i]*LenVal[i]/(fc+ft) + alpha1[i] =0 ; 

    subject to crushing_bot {i in lenline}:
    -x[3*i]*LenVal[i]/2 - x[3*i-1]*((fc-ft)*LenVal[i]/(2*(ft+fc)) -x[3*i-1]/(2*(fc+ft) )) - 0.5*fc*ft*LenVal[i]*LenVal[i]/(fc+ft) + alpha2[i] = 0; 
 
      maximize total_cost: 
        sum {i in cline} x[CInd[i]+1]*CVal[i] ;  
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
