% The hyper-Dirichlet distribution (Hankin 2010) with optimization (Dong and Tian 2013)
% Symbolic using mupad (may not work with maple)
classdef hdirichlet
%Design notes
% SetAccess
%  I try to balance between logic and automation. 
%  SetAccess will be public as long as there is no unacceptable breach of data consistency.
%  For those warned public SetAccess properties, the burden of maintaining data consistency is shared by the user
%=== INSTANCE MEMBERS Begin ===================================================================
    properties(SetAccess='private') % essential properties
        a; % a is the vector of ionic counts. Ensure: isnumeric(a) && all(a > -1)
        b; % b is the vector of unionic counts. Ensure: isnumeric(b)
        De;% De is a binary matrix, jth row indicates which components of p to be included in the summation to match the jth component of b; De is represented as of class 'logical'. Ensure: islogical(De) && all(size(De) == [length(b), length(a)])        
    end
    properties(SetAccess='private') % auxilliary symbolic properties for slicing algorithm
        x; % array of symbols, content like [x1, x2, ..., x20]
        tau; % array of symbols, content like [tau1,tau2,...,tau20]
        slices; % struct holding strings of polynomial equaltions
    end
    properties(SetAccess='public')  % setting & backpack
        setting = struct(... 
                      'TolX',1e-10, ...
                      'MaxIter', 1000, ...
                      'SaveIteration', true, ...
                      'ShowAlert', true ...
                   );
        backpack;
    end
    methods % THE CONSTRUCTOR
        function this = hdirichlet ... % constructing from only the vitals --> post: optimization ready
                 (ionic_counts, ...                % --> this.a
                  unionic_counts, ...              % --> this.b
                  unionic_event_pattern_matrix ... % --> this.De (length of column = length of ionic_counts)
                 )
            
            if nargin == 0
                this.a = [0;0];
                this.b = 0;
                this.De = false(2,1);
            elseif nargin == 1 % only ionic_counts ==> degenerate into Dirichlet
                this.a = ionic_counts(:);
                this.b = 0;
                this.De = repmat(logical(false),length(this.a),1); % doesn't matter true/false as b is 0
            elseif nargin == 2 % only unionic_counts & unionic_event_pattern_matrix
                disp('NB: Assuming the two inputs are for the unionic_counts and the unionic_event_pattern_matrix');
                this.De = logical(unionic_counts);
                this.b = ionic_counts(:);
                if size(this.De, 2) ~= length(this.b) % strict validation is delayed
                    this.De = this.De.';
                end
                this.a = zeros(size(unionic_event_pattern_matrix,1),1);
            elseif nargin == 3
                this.a = ionic_counts(:);   % force column vector
                this.b = unionic_counts(:); % force column vector
                this.De = logical(unionic_event_pattern_matrix); % force bit matrix
            end
            % Transpose this.De when obvious oriention error detected. Warning: When a and b have equal length, De's dimension error can't be detected
            if (length(this.a) ~= length(this.b)) ... % not square matrix
            && all(size(this.De) == [length(this.b),length(this.a)]) % and in correct dimensions but wrong position               
                this.De = this.De.';
                disp('NB: The composite indicator matrix, De, was input as size length(b) * length(a). I have assumed the intended input was the transpose.');                
            end
            
            assert(this.validate_vital(), ...
                   'Some of the vitals(a, b, and De) has FAILED validation');
            this = this.checkIsDirichlet();
            if size(this.De) * [1;-1] == 0
                disp('Caution: The composite indicator matrix, De, is a square matrix => Attention to its orientation.');
            end
            
        end
    end % THE CONSTRUCTOR
    methods % VALIDATOR
% VALIDATOR methods index
% validate_vital() as hdirichlet .............. wrapper of validators for a, b, De (the vitals)
% validate_vital_syntax() as logical........... halts execution flow whenever any (detectable) error in the syntax (data type, size, etc) of a,b,De (the vitals)
% validate_regularity() as logical............. check whether the current distribution is convergent on the boundaries. c.f. thm 26 (Criterion for uniform regularity)

        function res = validate_vital_syntax ... % validator for the types and sizes of vital properties
                (this ...
                )
            assert(isnumeric(this.a),   'a is not a numeric vector.');
            assert(all(this.a > -1),    'a is not all > -1');
            assert(isnumeric(this.b),   'b is not a numeric vector.');
            assert(islogical(this.De),  'De is not a logical matrix.');
            assert(all(size(this.De) == [length(this.a), length(this.b)]), ...
                                        'size(De) is not [length(a), length(b)]');
            res = true;
        end
        
        function res = validate_vital... % validator for vital properties
                 (this ...
                 )
            assert(this.validate_vital_syntax(), ...
                                        'Some vital parameters have wrong type or size.');
            res = true;
        end
        function res = validate_regularity_necessary ... % [~N*C^N] a quite tight necessary condition only
                 (this ...
                 )
             res = false;
             bigb = [this.a; this.b].';
             bigDe = [eye(length(this.a)), this.De];
             id = 1:(length(bigb));
             idnegb = id(bigb <= 0);
             for j = idnegb % individual order-1 frag
                 if sum(bigb(:,~any(bigDe(~bigDe(:,j),:),1))) < 0
                     disp(strcat('error @ ', num2str(j)));
                     return; % => must be irregular
                 end
             end
             
             res = true; % => probably regular
        end
        
        function res = validate_regularity_single_node ...
                 (this, ...
                  bitvec ...
                 )             
             bitvec = logical(bitvec(:));
             if all(bitvec) % => union is exhaustive => always regular
                 res = true;
                 return;
             end
             res = false;
             bigb = [this.a; this.b].';
             bigDe = [eye(length(this.a)), this.De];   
             assert(any(and(all(bigDe(bitvec,:)), ~any(bigDe(~bitvec,:)))));
             if sum(bigb(:,~any(bigDe(~bitvec,:),1))) < 0                 
                 return; % => must be irregular
             end
             res = true;
        end
        
        function res = validate_regularity_union_node ...
                 (this, ...
                  bitvecs ...
                 )             
             if size(bitvecs, 1) ~= length(this.a) && size(bitvecs, 2) == length(this.a)                 
                 bitvecs = bitvecs.';
             end             
             res = false;
             bitvec = any(bitvecs,2);
             if all(bitvec) % => union is exhaustive => always regular
                 res = true;
                 return;
             end
             [~,id] = ismember(bitvecs.',this.De.','rows');
             bigb = [this.a; this.b; sum(this.b(id))].';
             bigDe = [eye(length(this.a)), this.De, bitvec];
             if sum(bigb(:,~any(bigDe(~bitvec,:),1))) < 0                 
                 return; % => must be irregular
             end
             res = true;
        end
    end % VALIDATOR
    methods % OPTIMIZATION
% Optimization methods index
% weaver() as struct .............. The Weaver Algorithm. Returns a struct holding the eigenestimate and iteration path.
% greedyweaver() as struct ........ The Greedy Weaver Algorithm. Returns a struct holding the eigest and iteration path.
% mode() as struct ................ wrapper and hybrider of weaver() and greedyweaver(). Returns a struct holding eigest and iters.
% p() as array of doubles ......... wrapper and convenienter of mode(). Returns the eigenestimate vector.
% Tau() as array of doubles ....... Returns the local thicknesses of unionic slices
% Tau0() as double ................ Returns the local thickness of the ionic slice
% d() as array of doubles ......... Returns the ionic counts reconstruction deviations
% I() as positive integer ......... Returns the index of the component of d(x) with the largest ABSolute value
% sse as double ................... Returns sum of squared inoic counts reconstruction deviations
% guessInit() as array of double .. Returns guessed initial eigenestimate to start weaver() or greedyweaver()

        function ini = guessInit ... % Returns guessed initial eigenestimate to start weaver() or greedyweaver()
                 (this ...
                 )
            ini = this.a / sum(this.a);
        end
        function res = p ... % shortcut of mode, returns only the p
                 (this ...
                 )
             this.setting.SaveIteration = false;
             this.setting.ShowAlert = false;
             res = this.mode.p;
        end
        function res = mode ... % mode optimizaiton
                 (this, ...
                  method, ... % 'weaver','greedyweaver'
                  ini ... % starting point for iteration
                 )
            if nargin == 1
                method = 'try';
                ini = 'default';
            elseif nargin == 2                
                ini = 'default';
            else
                assert(strcmpi(method, 'try')||strcmpi(method, 'weaver')||strcmpi(method, 'greedyweaver'), ['method ', method, ' is not one of {''nr'', ''gm'', ''nm''}']);
            end
            
            if strcmp(method, 'try')
                try
                    if this.setting.ShowAlert
                        disp('trying weaver algorithm')
                    end
                    res = this.weaver(ini);
                    return;
                catch e_weaver
                    if this.setting.ShowAlert
                        disp('trying greedy weaver')
                    end
                    res = this.greedyWeaver(eval(e_weaver.message)); 
                    return;
                end
            end
            
            if strcmpi(method, 'weaver') % stub: weaver
                res = this.weaver(ini);
                return;                
            end
            
            if strcmpi(method, 'greedyweaver') % stub: reconstruction-conditional tying
                res = this.greedyWeaver(ini);
                return;                
            end
        
        end
        function res = Tau ... x-thickness of every unionic event
                 (this, ...
                  x ...
                 )
             res = this.b ./ (this.De.' * x);
        end        
        function res = Tau0 ... x-thickness of the ionic slice
                 (this, ...
                  x ...
                 )
             res = (x.' * (this.a - diag(x) * (~this.De) * this.Tau(x))) / (x.' * x);
        end                
        function res = d ... (ionic counts reconstructed at x) - the actual ionic counts
                 (this, ...
                  x ...
                 )
             res = x .* (this.Tau0(x) + (~this.De) * this.Tau(x)) - this.a;
        end        
        function res = sse ...
                 (this, ...
                  x ...
                 )
             res = (this.d(x)).' * this.d(x);
        end        
        function res = I ... returns the index of the component of d(x) with the largest absolute value
                 (this, ...
                  x ...
                 )
             dd = this.d(x);
             [~, res] = max(abs(dd(1:length(dd)-1)));
        end        
        function res = weaver ...
                 (this, ...
                  ini ...
                 )
            if nargin < 2
                ini = this.guessInit;
            end
            if abs(sum(ini) - 1) > eps * 100 || ~all(ini > 0)        
                if this.setting.ShowAlert
                    disp('User''s initial value for Weaver has a problem --> I''ll use default initial value to start iterations!');                
                end
                ini = this.a / sum(this.a);
            end
            x = ini; 
            smallestSSE = this.sse(x);
            bestx = x;
            iter = x;
            itercount = 0;
            while this.sse(x) > this.setting.TolX ^ 2 % squared because of `squared' sum of errors
                x = this.a ./ ((~this.De) * this.Tau(x) + this.Tau0(x));
                if any(isnan(x)) || itercount > this.setting.MaxIter || max(x) == Inf 
                    e_weaver = MException('MATLAB:hdirichletMaxIterReached',mat2str(bestx));
                    if this.setting.ShowAlert
                        warning('Weaver failed: Either this.setting.MaxIter is reached or divergence has occured. The best estimate is captured in the message. You may want to validate regularity.');
                    end
                    throw(e_weaver); % the e_weaver.message string captures the best result before failure
                end
                x = x / sum(x);
                if this.setting.SaveIteration
                    iter = [iter, x];
                end
                if this.sse(x) < smallestSSE
                    smallestSSE = this.sse(x);
                    bestx = x;
                end
                itercount = itercount + 1;                
            end
            res = struct('p',bestx,'Tol_sse', this.setting.TolX ^ 2, 'lnUdensity', this.lnUdensity(bestx),'iter',iter.');            
        end        
        function res = greedyWeaver ...
                 (this, ...
                  ini ...
                 )
            if nargin < 2
                ini = this.guessInit;
            end
            if abs(sum(ini) - 1) > eps * 100 || ~all(ini > 0)       
                if this.setting.ShowAlert
                    disp('User''s initial value for Greedy Weaver has a problem --> I''ll use default initial value to start iterations!');                
                end
                ini = this.a / sum(this.a);
            end
            x = ini;
            n = length(this.a);  
            smallestSSE = this.sse(x);
            bestx = x;
            iter = x;
            itercount = 0;
            while this.sse(x) > this.setting.TolX ^ 2
                i = this.I(x);
                u = [0,x(i), 1.05 * x(i)];
                temp = x;
                temp(i) = u(3);
                temp(n) = 0; 
                temp(n) = 1 - sum(temp);
                dp = this.d(x);
                dtemp = this.d(temp);
                v = [-this.a(i), dp(i), dtemp(i)];
                ga = v(1);
                al = (u(3) * (v(2) - v(1)) - u(2) * (v(3) - v(1))) / (u(2) * u(2) * u(3) - u(2) * u(3) * u(3));
                be = (-u(3) * u(3) * (v(2) - v(1)) + u(2) * u(2) * (v(3) - v(1))) / (u(2) * u(2) * u(3) - u(2) * u(3) * u(3));
                x(i) = (-be + sqrt(be * be - 4.0 * al * ga)) / (2.0 * al);
                x(n) = 0;
                x(n) = 1 - sum(x);
                if this.setting.SaveIteration
                    iter = [iter, x];
                end
                if this.sse(x) < smallestSSE
                    smallestSSE = this.sse(x);
                    bestx = x;
                end
                itercount = itercount + 1;
                if any(isnan(x)) || itercount > this.setting.MaxIter || max(x) == Inf
                    e_greedyWeaver = MException('MATLAB:hdirichletMaxIterReached',mat2str(bestx));
                    if this.setting.ShowAlert
                        warning('Greedy Weaver failed: Either this.setting.MaxIter is reached or divergence has occured. The best estimate is captured in the message. You may want to validate regularity.');
                    end
                    throw(e_greedyWeaver); % the e_greedyWeaver.message string captures the best result before failure
                end
            end
            res = struct('p',bestx,'Tol_sse', this.setting.TolX ^ 2, 'lnUdensity', this.lnUdensity(bestx),'iter',iter.');
        end
    end % OPTIMIZATION
    methods % SYMBOLICS
%SYMBOLICS methods
% constructVitalSymbols() as hyperdirihclet ....... construct the vital symbols
% udensitySym() as expression ..................... returns symbolic expression of the un-normalized density
% generateSymsDeclarationString() as string ....... generate an evaluation-ready expression for adding/clearing the symbols to/from the global scope        

        function this = constructVitalSymbols ...
                (this ...
                )
            assert(~isempty(this.a));
            for i = 1:length(this.a)
                this.x = [this.x; sym(['x',num2str(i)],'real')];
            end
            this.tau = sym('tau0');
            if ~this.isDirichlet
                for i = 1:length(this.b)
                    this.tau = [this.tau; sym(['tau',num2str(i)])];
                end
            end
        end
        function res = udensitySym ...
                (this ...
                )
            assert(~isempty(this.a));
            assert(~isempty(this.x));
            res = prod(this.x .^ this.a);
            if ~this.isDirichlet
                res = res * prod((this.De' * this.x).^ this.b);
            end
        end
        
        function res = generateSymsDeclarationString ...
                (this, ...
                 command ... % usually either 'syms' or 'clear'
                )
        % In global scope execute: eval(myHdirInstance.generateSymsDeclarationString)            
            assert(~isempty(this.a));
            if nargin == 1
                command = 'syms';
            end
            res = command;
            for i = 1:length(this.x)
                res = [res, ' x', num2str(i)];
            end
            res = [res, ' real; '];
            res = [res, command];
            for i = 1:length(this.tau)
                res = [res, ' tau', num2str(i-1)];
            end
            res = [res, ';'];
                
        end
        function [this, sol] = sliceIntoEquations ...
                (this ...                 
                )
            n = length(this.a);
            this.slices = sum(this.x) - 1;                             
            if this.isDirichlet
                for i = 1:n
                    eq = this.x(i) * this.tau(1) - this.a(i);
                    this.slices = [this.slices, eq];
                end
                return;
            end
            
            m = length(this.b);
            im = 1:m;            
            for i = 1:n % part I: atoms: n equations each of shape x1 * (tau0 + tau2 + tau4)
%                 if any(im(s==i))
                    eq = this.x(i) * (this.tau(1) + sum(this.tau(im(~this.De(i,:))+1))) - this.a(i);
%                 else
%                     eq = this.x(i) * this.tau(1);
%                 end
                this.slices = [this.slices, eq];
            end
            this.slices = [this.slices, ((this.De.' * this.x) .* this.tau(2:length(this.tau)) - this.b).'];            
            
            if nargout == 2
                sol = solve(this.slices);
            end
        end
        function [x, tau] = solveSlices ...
                (this ...
                )
            m = length(this.b);
            if ~this.isDirichlet
                m = m + 1;
            end
            assert(length(this.slices) == length(this.a) + m) % considering adding a property to indicate whether the slices are properly prepared
            
            sol = solve(this.slices);
            x = [];
            for i = 1:(length(this.a))
                eval(['x=[x,sol.x', num2str(i), '];']);
            end
            if nargout == 2                
                tau = sol.tau0;
                if m >= 2
                    for i = 1:(m-1)
                        eval(['tau =[tau,sol.tau', num2str(i),'];']);
                    end
                end
            end
        end
    end % SYMBOLICS
    methods % UTILITY 
% Utility methods index
% isDirichlet() as logical ....................... return true if this hyper-dirichlet object is a dirichlet one
% checkIsDirichlet() as hdirichlet ............... a macro that wraps isDirichlet() and do the trimming on (b, De)
% reduce() as hdirichlet ............... collect the powers for same composites and trim all-zero and all-one columns and order the delta matrix
        function res = isDirichlet ... % post: return true if this hyper-dirichlet object is a dirichlet one
                 (this ...
                 )
            res = false;
            if all(this.b == 0) 
                res = true;
                return;
            end
            
            if all(all(this.De == 1)) || all(all(this.De == 0))
                res = true;
                return;
            end
            
            if all(this.De(repmat((this.b ~= 0), size(this.De,1), 1)) == 1) || ...
               all(this.De(repmat((this.b ~= 0), size(this.De,1), 1)) == 0)
               res = true;
               return;
            end
        end
        
        function this = checkIsDirichlet ... % wrapper for Dirichlet degeneracy checking
                 (this ...
                 )
            if this.isDirichlet()
%                 res = input('This is a dirichlet distribution, degenerate the composite powers(b) and composite indicator matrix(De)? [y/n]:','s');
% (2012-07-05) jdong removed the preceding line and use the following line, for better automation flow
                res = 'y'; disp('This is a dirichlet distribution. I will degenerate the composite powers(b) and composite indicator matrix(De).');
                
                if (res == 'y') || (res == 'Y')
                    this.b = 0;
                    this.De = true(length(this.a),1); % use true for later compatibility
                end
            end
        end
        
        function this = reduce ... % [~N] look for same columns in De and combine them, also triming all-zero and all-one, and order the Delta matrix
                 (this ...
                 )
             Deb = sortrows([this.De.',this.b,(1:length(this.b)).']);             
             if all(Deb(1,1:length(this.a)) == 0) % to trim heading allzero rows in De
                 i = 1;
                 while all(Deb(i,1:length(this.a)) == 0)
                     i = i + 1;
                 end
                 Deb = Deb(i:end,:);
             end
             if all(Deb(end,1:length(this.a)) == 1) % to trim tailing allone rows in De
                 i = size(Deb,1);
                 while all(Deb(i,1:length(this.a)) == 1)
                     i = i - 1;
                 end
                 Deb = Deb(1:i,:);
             end            
             Deb = Deb(Deb(:,length(this.a)+1) ~= 0,:);  % remove rows where b = 0
             [n,q] = size(Deb);
             newb = zeros(n,1);
             newDe = zeros(n,q-2);
             newb(1) = Deb(1,q-1);
             newDe(1,:) = Deb(1,1:(q-2));
             j = 1;
             if n > 1
                 for i = 2:n
                     if any(Deb(i,1:(q-2)) ~= newDe(j,:))
                         j = j + 1;
                         newDe(j,:) = Deb(i,1:(q-2));
                     end
                     newb(j) = newb(j) + Deb(i,q-1);
                 end
             end
             this = hdirichlet(this.a, newb(1:j), newDe(1:j,:).');
        end
    end % UTILITY
    methods % GRAPHICS 
% Graphics methods index
%plotStrands() as struct ................ plots length(a) + length(b) number of log strands in the summation of log-density; also allowing user specified thetas to be plotted
        function res = plotStrands ... %plots length(a) + length(b) number of log strands in the summation of log-density; also allowing user specified thetas to be plotted
                 (this, ...                  
                  p0, ... % accepts a matrix, 0, or a positive integer(randomly generate this number of p0)
                  xlength ...
                 )
            if nargin == 1
                p0 = ones(1, length(this.a)) ./ length(this.a);
            end
            if nargin >= 2
                if ~all(size(p0) == 1)
                    assert((size(p0,1) == length(this.a)) || (size(p0,2) == length(this.a)), 'p0 has wrong size');                
                    assert(all(p0(:) >= 0) && (all(abs(sum(p0) - 1) < this.setting.TolX) || all(abs(sum(p0.') - 1)  < this.setting.TolX)), 'p0 is not on the simplex');
                    if (...
                        (size(p0, 1) == length(this.a) && all(abs(sum(p0) - 1) < this.setting.TolX)) ...
                    && ~(size(p0, 2) == length(this.a) && all(abs(sum(p0.') - 1)  < this.setting.TolX)))
                        p0 = p0.';                    
                    end
                else
                    assert(p0 >= 0);
                    p0 = hdirichlet.simplexUnifRnd(length(this.a), p0); %after this line p0 is a matrix
                end
                    
            end
            if nargin < 3
                xlength = 50;
            end
            
            lenp0 = size(p0,1);
            x0 = [p0, p0 * this.De];
            y0 = log(x0) .* repmat([this.a; this.b].', lenp0, 1);
            res = struct('x0', x0, 'y0', y0, 'p0', p0);  
            
            hold on;
            x = linspace(0,1,xlength);
            y = [this.a; this.b] * log(x);
            plot(x,y,':');
            
%             mfc = unifrnd(1, 3*lenp0, lenp0, 3);
%             marks = 'sod^vph'; % marker shapes
%             mecs = 'rmgycbk'; % marker edge colors
%             for i=1:lenp0
%                 mfcc = mfc(i,:) / sum(mfc(i,:));
%                 mfcc = max(mfcc - median(mfcc),0) * 5 + median(mfcc);
%                 mfcc = mfcc / sum(mfcc);                
%                 plot(x0(i,:),y0(i,:),marks(floor(unifrnd(1,length(marks)-0.01))),'LineWidth', 2, 'MarkerSize',floor(unifrnd(5,8.99)), 'MarkerEdgeColor',mecs(floor(unifrnd(1,length(mecs)-0.01))), 'MarkerFaceColor', mfcc);
%                 
%             end
            plot(x0',y0','s','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','g');
            legend('show');            
            hold off;
        end
        
    end % GRAPHICS
    methods % DENSITY
% Density methods index
% lnUdensity as double ............... ln (base e) un-normalized density 
% lgUdensity as double ............... lg (base 2) un-normalized density 
% udensity as double ................. un-normalized density
% densityIntegration as double ....... normalizing constant = integration of the un-normalized density
% lgDensityIntegration as double ..... integration of the log-un-normalized density
        function res = lnUdensity ...
                 (this, ...
                  p ...
                 )
            if size(p,2) == length(this.a) ...
            && size(p,1) ~= length(this.a)
                p = p.';
            end
            id = (this.a ~= 0);
            res = this.a(id).' * log(p(id,:)) + this.b.' * log(this.De.' * p);
            res = res(:);
        end
        
        function res = lgUdensity ...
                 (this, ...
                  p ...
                 )
            res = this.lnUdensity(p) ./ log(2);
        end
        
        function res = udensity ...
                 (this, ...
                  p ...
                 )
            res = exp(this.lnUdensity(p));
        end
                
        function res = lgDensityIntegration ... % Caution: Expecation(uniform)[lg(Unnormalized density)] < lg(Expectation(uniform)[Unormalized density])
                (this, ...
                 sizePerMonteCarlo, ...
                 nMonteCarlo ...
                )
            if nargin == 1
                sizePerMonteCarlo = 5000;
                nMonteCarlo = 100;
            elseif nargin ==2
                nMonteCarlo = 100;
            end
            
            res = 0;
            for i = 1:nMonteCarlo
                res = res + sum(this.lgUdensity(hdirichlet.simplexUnifRnd(sizePerMonteCarlo,length(this.a)))) / sizePerMonteCarlo; 
            end
            res = res / nMonteCarlo;
        end  
        
        function res = densityIntegration ... % !!!dangerous version!!! do not use
                (this, ...
                 sizePerMonteCarlo, ...
                 nMonteCarlo ...
                )
            warning('!!! densityIntegration is currently dangerous, result often very wrong !!! do not use');
            if nargin == 1
                sizePerMonteCarlo = 5000;
                nMonteCarlo = 100;
            elseif nargin ==2
                nMonteCarlo = 100;
            end
            
            res = 0;
            for i = 1:nMonteCarlo
                res = res + sum(this.udensity(hdirichlet.simplexUnifRnd(sizePerMonteCarlo,length(this.a)))) / sizePerMonteCarlo; 
            end
            res = res / nMonteCarlo;
        end 
    end % DENSITY
    methods(Static) %Class Methods
% Static methods index
% parseVitalFromPowerSet as arraylist ...... converting a vector of powers, ordered to match the binary expansion of indicators, from 00..01 to 11..10, total 2^n powers, into the 3 vitals that is ready to construct an instance of hdirichlet
% embedVitalIntoPowerset as double array ... inverse function of parseVitalFromPowerSet
% simplexUnifRnd as double array ........... returns a uniformly distributed point on the simplex
% powerset as logical matrix ............... returns 2^n rows of bits indicating subsets of a finite set of n elements
% parseVitalFromWeavingGrid ................ converting the weavingGrid in the format of Fig 2.3 of Dong and Tian (2013) excluding the ?s into a, b, and De, the three essential parameter of an hdirichlet object.
% embedVitalIntoWeavingGrid ................ inverse function of parseVitalFromWeavingGrid

        function [a,b,De] = parseVitalFromPowerset ... % converting a vector of powers, ordered to match the binary expansion of indicators, from 00..00 to 11..11, total 2^n powers, into the 3 vitals that is ready to construct an instance of hdirichlet
                 (powers ...
                 )
            % powers are assumed in the order of binary expansion starting 00..00, ending 11..11
            % eg powers = [1,2,3,4] corresponds to 2-digit binary bits [00,01,10,11]
            disp('I assume that you are pasting from the R package hyperdirichlet format, authored by Hankin (2010) where the ionic counts are decremented by 1. ---> Effectively this method is parameter converter between the R::hyperdirichlet and my Matlab hdirichlet class.');
            assert(isnumeric(powers) && any(size(powers) == 1), 'input "powers" is not a numeric vector.');            
            n = length(powers);
            assert( (n >= 2) && (2^(log2(n)) == n), 'The length of the input "powers" is not 2^n where n >= 2.');
            powers = powers(:);
            powers = powers(2:n-1);% trim the useless first and last elem of powers
            n = n - 2;
            na = log2(n+2); %length of a
            ia = bitshift(1,(0:(na-1))); %int32 position of ionic_counts in powers
            a = flipud(powers(ia)) - 1;
            ib = true(1,n);
            ib(ia) = false; %logical subscript vector of unionic_counts in powers
            b = powers(ib);
            t = 1:n; 
            De = (char(dec2bin(t(ib))) == '1').';
        end
        function powers = embedVitalIntoPowerset ... % inverse function of parseVitalFromPowerSet
                 (a, ...
                  b, ...
                  De ...
                 )
             disp('This works as parameter converter from my Matlab hdirichlet class to Hankin(2010) R::hyperdirichlet class ---> the ionic counts there is incremented by 1.')
             if (~all(size(De) == [length(a), length(b)])) ...
             && (all(size(De) == [length(b), length(a)]))
                De = De.';
             end
             a = flipud(a(:))+1;
             b = [flipud(a(:)); b(:)];
             
             De = [eye(length(a)), De];
             t1 = num2str(De).';
             t2 = t1(t1(:,1)~=' ',:);
             t3 = bin2dec(t2);
             assert(length(t3) == length(b));            
             powers = zeros(2^length(a) - 2, 1);
             powers(t3) = b;
             powers = [0;powers;0];
        end
        function res=powerset(n)
            res = false(2^n,n);
            for i=1:n
                res(:,i) = repmat([false(2^(i-1),1);true(2^(i-1),1)],2^(n-i), 1);        
            end
        end
        function res = simplexUnifRnd ... % a row is a uniformly distributed point on the simplex
                 (sampleSize, ...
                  length ... %dim = length-1
                 )
             % fact: difference sequence of [0, unif[0,1], 1] is unif on simplex
             assert(length >=2 && sampleSize >= 1);
             res = diff([zeros(1,sampleSize);sort(unifrnd(0,1,length-1,sampleSize));ones(1,sampleSize)]).';
        end
        function [a,b,De] = parseVitalFromWeavingGrid ...
                 (weavingGrid ... % automatically check whether top row is all zero
                 )
             disp('If you are pasting the weaving grid from my website, pls fill the top-right and bottom-right grid with an arbitrary number.');
             dim = size(weavingGrid);             
             a = weavingGrid(dim(1),1:(dim(2)-1)).';
             firstRow = 1;
             if max(weavingGrid(1,1:dim(2))) == 0 && min(weavingGrid(1,1:dim(2))) == 0
                 firstRow = 2;
             end
             b = weavingGrid(firstRow:(dim(1)-1), dim(2));
             De = weavingGrid(firstRow:(dim(1)-1), 1:(dim(2)-1)).';
        end
        function weavingGrid = embedVitalIntoWeavingGrid ... % with first row all zero
                 (a, ...
                  b, ...
                  De ...
                 )
             weavingGrid = [zeros(1,length(a)+1);[De.',b;a.',0]];
        end
    end %Class Methods
end
