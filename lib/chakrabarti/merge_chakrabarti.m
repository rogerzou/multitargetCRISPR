A = readtable("chakrabarti_mm2.xlsx");
B = readtable("chakrabarti_mm3.xlsx");
[~,ia,ib] = intersect(A(:,1),B(:,1));
C = [A(ia,:),B(ib,2:end)];
writetable(C, "chakrabarti_mm.csv")
