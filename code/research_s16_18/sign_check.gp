\\ Sign of PARI's plus modular symbol at the cusp 0 versus L(E,1)/omega_1, on rank-0 curves.
\\ Mazur-Tate normalize [0]^+ = +L(E,1)/Omega^+; PARI's msfromell gives the opposite sign (ratio -1).
\\ usage: gp -q sign_check.gp   (or read() from cypari2)
t(a)={my(E=ellinit(a),Ms=msfromell(E,1),L1=lfun(lfuncreate(E),1),s=mseval(Ms[1],Ms[2],[0,oo]),w=real(E.omega[1]));
 printf("%s disc%s L(1)=%.4f  [0,oo]_PARI=%s  L/omega1=%.4f  ratio=%.3f\n",a,if(sign(E.disc)>0,">0","<0"),L1,s,L1/w,s/(L1/w))}
t([0,-1,1,-10,-20]); t([1,0,1,4,-6]); t([1,-1,1,-1,-14]); t([0,1,1,-1,0]); t([0,1,1,-23,-50]); t([1,1,1,-10,-10]); t([1,-1,0,-1,0]);
