function [final_Rscore]=GRMDA(interaction,nd,nm,sd,sm)

[ua,ba,va]=svds(interaction,180);
Ar=ua*sqrt(ba);
Ad=(va)*(sqrt(ba))';

[um,bm,vm]=svds(sm,220);
Fr=um*sqrt(bm);
[ud,bd,vd]=svds(sd,170);
Fd=ud*sqrt(bd);


Br=PLS(Fr,Ar);
Bd=PLS(Fd,Ad);
final_Rscore=Fr*Br*Bd'*Fd';
end
%D=1./(1+exp(-C));

%allresult(miRNA,disease,A,C);