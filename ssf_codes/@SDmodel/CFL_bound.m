function tCFL = CFL_bound(obj,dz,qhat,temperature,vbmax,phiclog)
            arguments
                obj (1,1) SDmodel
                dz = 1e-3;
                qhat = 7.2/0.4;
                temperature = 273;
                vbmax = 1e3;
                phiclog = .8;
            end
            [w_v, w_a, w_b, w_s] = deal(zeros(5,1));

            vfmax = (qhat + vbmax*phiclog)/(1 - phiclog);

            rhoP = max(obj.particles.density);
            alfaP = max(obj.particles.dispersivity);
            alfaL = max(obj.particles.dispersivity);

            att_max = max(obj.particles.attachment_rate);
            trP_max = max(obj.particles.transfer_rate);
            trL_max = max(obj.particles.transfer_rate);

            SP = obj.ecological.stoichiometric_matrix{obj.particles.Name,:};
            SL = obj.ecological.stoichiometric_matrix{obj.liquids.Name,:};
            SL(SL > 0) = 0;
            K = obj.ecological.kinetic_parameters(obj.liquids.Name,:);
            mu = obj.ecological.kinetic_parameters{"Reaction Rate",:}.*...
                (obj.ecological.kinetic_parameters{"Temperature Correcting Factor",:}.^(293 - temperature));
              

            w_v([1 2 3]) = 2*vbmax;
            w_v([4 5]) = 2*vfmax;

            w_a([1 2 3]) = 0;
            w_a(4) = 2*vfmax*alfaP*(1 + 1/(1 - phiclog));
            w_a(5) = 2*vfmax*alfaL*(1 + 1/(1 - phiclog));

            w_b(1) = obj.detachment(vfmax);
            w_b(2) = att_max + trP_max/obj.beta;
            w_b(3) = max(trL_max/obj.beta,1/obj.tau);
            w_b(4) = att_max*(1 + phiclog) + trP_max/obj.beta*(phiclog/(1 - phiclog));
            w_b(5) = trL_max/obj.beta*(phiclog/(1 - phiclog));

            w_s([1 2 3]) = max([SP(1,3)*mu(3), SP(2,4)*mu(4), SP(3,5)*mu(5)/K{5,5}]);
            w_s(4) = 2*max(sum(abs(SL)*(rhoP*(mu(:).')./K{:,:}),2));
            w_s(5) = max(sum(abs(SL)*(rhoP*(mu(:).')./K{:,:}),2));

            w = w_v/dz + w_a/dz^2 + w_b + w_s;

            tCFL = min(1./w);
        end