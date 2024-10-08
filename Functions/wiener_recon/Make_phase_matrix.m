function [phase_matrix] = Make_phase_matrix(initial_phase,spacing,regul)
    deph1 = initial_phase(1,1);    
    deph2 = initial_phase(2,1);
    deph3 = initial_phase(3,1);
    deph4 = initial_phase(1,2);    
    deph5 = initial_phase(2,2);
    deph6 = initial_phase(3,2);

    phase_matrix= [1 1 1 1 1;
        exp(2*1i*regul*(spacing(1)/sum(spacing)))  exp(1i*regul*(spacing(1)/sum(spacing))) 1  exp(-1i*regul*(spacing(1)/sum(spacing))) exp(-2*1i*regul*(spacing(1)/sum(spacing)));
        exp(2*1i*regul*((spacing(1)+spacing(2))/sum(spacing))) exp(1i*regul*((spacing(1)+spacing(2))/sum(spacing))) 1  exp(-1i*regul*((spacing(1)+spacing(2))/sum(spacing))) exp(-2i*regul*((spacing(1)+spacing(2))/sum(spacing)));
        exp(2*1i*regul*((spacing(1)+spacing(2)+spacing(3))/sum(spacing))) exp(1i*regul*((spacing(1)+spacing(2)+spacing(3))/sum(spacing))) 1  exp(-1i*regul*((spacing(1)+spacing(2)+spacing(3))/sum(spacing))) exp(-2i*regul*((spacing(1)+spacing(2)+spacing(3))/sum(spacing)));
        exp(2*1i*regul*((spacing(1)+spacing(2)+spacing(3)+spacing(4))/sum(spacing))) exp(1i*regul*((spacing(1)+spacing(2)+spacing(3)+spacing(4))/sum(spacing))) 1  exp(-1i*regul*((spacing(1)+spacing(2)+spacing(3)+spacing(4))/sum(spacing))) exp(-2i*regul*((spacing(1)+spacing(2)+spacing(3)+spacing(4))/sum(spacing)));];
    phase_matrix1 = phase_matrix.*repmat([exp(1i*deph1),exp(1i*deph4),1,exp(-1i*deph4),exp(-1i*deph1)],[size(phase_matrix,1),1,1]);
    phase_matrix2 = phase_matrix.*repmat([exp(1i*deph2),exp(1i*deph5),1,exp(-1i*deph5),exp(-1i*deph2)],[size(phase_matrix,1),1,1]);
    phase_matrix3 = phase_matrix.*repmat([exp(1i*deph3),exp(1i*deph6),1,exp(-1i*deph6),exp(-1i*deph3)],[size(phase_matrix,1),1,1]);
    phase_matrix(:,:,1) = inv(phase_matrix1);
    phase_matrix(:,:,2) = inv(phase_matrix2);
    phase_matrix(:,:,3) = inv(phase_matrix3);
    
end