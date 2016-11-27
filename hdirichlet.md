<pre style="color:#FFF;font-size:110%;">
% 1: How to input data

>> a=[23,41,40,17].', b=[12,8,24,14,-22,-5,-3].', De=[1,1,0,0;0,0,1,1;1,0,1,0;0,1,0,1;0,1,1,0;1,0,0,1;1,0,1,1].'
a =
    23
    41
    40
    17
b =
    12
     8
    24
    14
   -22
    -5
    -3
De =
     1     0     1     0     0     1     1
     1     0     0     1     1     0     0
     0     1     1     0     1     0     1
     0     1     0     1     0     1     1
     
% Alternatively, if you would like to paste from the above Excel Web app, 
%   you may use the class method hdirichlet.parseVitalFormWeavingGrid that takes the weaving grid matrix as input.
% You must fill the top-right and bottom-right cell of the weaving grid with arbitrary numbers to make it a matrix
>> [a,b,De]=hdirichlet.parseVitalFromWeavingGrid([0	0	0	0	0
1	1	0	0	12
0	0	1	1	8
1	0	1	0	24
0	1	0	1	14
0	1	1	0	-22
1	0	0	1	-5
1	0	1	1	-3
23	41	40	17	0])
If you are pasting the weaving grid from my website, pls fill the top-right and bottom-right grid with an arbitrary number.
a =
    23
    41
    40
    17
b =
    12
     8
    24
    14
   -22
    -5
    -3
De =
     1     0     1     0     0     1     1
     1     0     0     1     1     0     0
     0     1     1     0     1     0     1
     0     1     0     1     0     1     1

% Finally, if you would like to paste from R::hyperdirichlet of Hankin(2010), you may use the class method hdirichlet.parseVitalFromPowerset that takes the powerset vector as input.
% Note that the first and the last components of the powerset vector are useless therefore can be arbitrary.
>> [a,b,De]=hdirichlet.parseVitalFromPowerset([ 0    18    41     8    42    14   -22     0    24    -5    24    -3    12     0     0     0])
I assume that you are pasting from the R package hyperdirichlet format, authored by Hankin (2010) where the ionic counts are decremented by 1. ---> Effectively this method is parameter converter between the R::hyperdirichlet and my Matlab hdirichlet class.
a =
    23
    41
    40
    17
b =
     8
    14
   -22
     0  &lt;--- not affecting result, will be trimmed by "reduce"
    -5
    24
    -3
    12
     0  &lt;--- 
     0  &lt;--- 
De =
     0     0     0     0     1     1     1     1     1     1
     0     1     1     1     0     0     0     1     1     1
     1     0     1     1     0     1     1     0     0     1
     1     1     0     1     1     0     1     0     1     0
>> h1=hdirichlet(a,b,De).reduce
h1 = 
  hdirichlet

  Properties:
           a: [4x1 double]
           b: [7x1 double]
          De: [4x7 logical]
           x: []
         tau: []
      slices: []
     setting: [1x1 struct]
    backpack: []
  Methods

% If you want to convert to the powerset format that R::hyperdirichlet reads, you may want to use the class method hdirichlet.embedVitalIntoPowerset
>> powerset = hdirichlet.embedVitalIntoPowerset(h1.a,h1.b,h1.De)
This works as parameter converter from my Matlab hdirichlet class to Hankin(2010) R::hyperdirichlet class ---> the ionic counts there is incremented by 1.
powerset =
     0
    18
    41
     8
    42
    14
   -22
     0
    24
    -5
    24
    -3
    12
     0
     0
     0




% 2: find the mode of the hyperdirichlet distribution

>> h1.p                            % request for the mode of the hyperdirichlet without any detail
ans =
      0.22879
       0.3126
      0.31534
      0.14328
>> h1.udensity(h1.p)               % unnormalized hyperdirichlet density
ans =
   1.6704e-80
>> h1.lnUdensity(h1.p)
ans =
      -183.69
>> res=h1.mode                     % request for the mode of the hyperdirichlet with details
trying weaver algorithm
User's initial value for Weaver has a problem --> I'll use default initial value to start iterations!
res = 
             p: [4x1 double]
       Tol_sse: 1e-20
    lnUdensity: -183.69            % ln(unnormalized hyperdirichlet density)
          iter: [22x4 double]
>> res.iter                        % iteration history
ans =
      0.19008      0.33884      0.33058       0.1405
      0.22897      0.31222      0.32129      0.13751
      0.22739      0.31427      0.31481      0.14354
      0.22892      0.31229      0.31576      0.14302
      0.22871      0.31273      0.31524      0.14332
       0.2288      0.31256      0.31537      0.14326
      0.22878      0.31261      0.31533      0.14328
      0.22879      0.31259      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879      0.31259      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
      0.22879       0.3126      0.31534      0.14328
>> h1.lnUdensity(res.iter)         % natural logarithm of the unormalized density
ans =
      -184.22
      -183.72
       -183.7
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69
      -183.69

>> h1.validate_regularity_necessary
ans =
     1                             % 1=regular (0=irregular): my weak method of validation does not indicate boundary irregularity here




% 3: symbolic features:

>> h1=h1.constructVitalSymbols     % initialize symbols for x's and tau's
h1 = 
  hdirichlet

  Properties:
           a: [4x1 double]
           b: [7x1 double]
          De: [4x7 logical]
           x: [4x1 sym]
         tau: [8x1 sym]
      slices: []
     setting: [1x1 struct]
    backpack: []
  Methods
>> h1.udensitySym                  % unnormalized hyper-dirichlet density in symbols
ans =
(x1^23*x2^41*x3^40*x4^17*(x1 + x2)^12*(x1 + x3)^24*(x2 + x4)^14*(x3 + x4)^8)/((x1 + x4)^5*(x2 + x3)^22*(x1 + x3 + x4)^3)
>> h1=h1.sliceIntoEquations        % the trivial slicing algorithm
h1 = 
  hdirichlet

  Properties:
           a: [4x1 double]
           b: [7x1 double]
          De: [4x7 logical]
           x: [4x1 sym]
         tau: [8x1 sym]
      slices: [1x12 sym]
     setting: [1x1 struct]
    backpack: []
  Methods
>> h1.slices.'                      % take a look at the TSA polynomials
ans =
                      x1 + x2 + x3 + x4 - 1
        x1*(tau0 + tau1 + tau2 + tau3) - 23
 x2*(tau0 + tau1 + tau4 + tau5 + tau6) - 41
        x3*(tau0 + tau2 + tau4 + tau7) - 40
        x4*(tau0 + tau3 + tau5 + tau7) - 17
                         tau1*(x3 + x4) - 8
                        tau2*(x2 + x4) - 14
                        tau3*(x2 + x3) + 22
                         tau4*(x1 + x4) + 5
                        tau5*(x1 + x3) - 24
                    tau6*(x1 + x3 + x4) + 3
                        tau7*(x1 + x2) - 12
>> h1.solveSlices                   % use commutative algebra to compute the root: total 19 roots found!! but only 1 on the simplex
ans =
[ -0.33297519240741753923340901805548, 0.50555334212384766055509386914727, 0.58092862671303900935018814828162,  0.24649322357053086932812700062659]
[ -0.75490704896460951373738713146565, 0.86049263351206190850781288829984, 0.17828467641654015224176766150622,  0.71612973903600745298780658165953]
[   0.5325926653846803367384146845042,  1.0221895698739220481601204821859, -1.2585944892941262406653503590478,  0.70381225403552385576681519235766]
[  0.22878819849569803604866673516895, 0.31259515443086074305433793267113, 0.31533654500403858181569688921201,  0.14328010206940263908129844294791] &lt;---
[  0.29537362604644408072857603339062,  0.5954489671981731333234854221725, 0.54186871468681732478232862234131, -0.43269130793143453883439007790443]
[   1.1717960849441753325859346224729, -1.0738633147493503048993021879178, 0.79579012315399958413446486160161,  0.10627710665117538817890270384323]
[   2.3963118410683933238969956938074,  1.0275786343361183328881049514765, -1.6039480284534258244027845333638, -0.81994244695108583238231611191907]
[  0.64860731534083353082475877677972, 0.82597826913866226253889442266754, 0.22653358371202584657566581496306, -0.70111916819152163993931901441014]
[  0.16355256245391517359785254903412, -1.0719636063448895171697333967193, 0.73314603457254117969495173806529,   1.1752650093184331638769291100357]
[  0.23271919799671427318666918172126, 0.66945706276035377561103393402126, 0.61442681658701841246932258647283, -0.51660307734408646126702570221535]
[ -0.62464031905239847400869040498613, 0.22472403969435838131162539135242, 0.82524622501246838292160032127386,  0.57467005434557170977546469235984]
[  0.12669158379508675413347843833616, 0.81911744133779593310777071267012,  -1.097059095012158700517112187182,   1.1512500698792760132758630363845]
[ -0.93131298787253830060388948756217,  1.0345545668573456315060811543495, 0.36018124278370698229209789119648,   0.5365771782314856868057104420162]
[ -0.80486183552268406215555451215111,  1.0328645972723320762034711017554,  1.5896838559143440104178708458078, -0.81768661766399202446578743541141]
[  0.47107580289772266381375265951619,  1.0552444196706162126068531519256, 0.40957547359347258760938556193538, -0.93589569616181146402999137338324]
[   1.2109666774583491136877669233189,  0.7166432817592774506017465007286, -1.0372488610004935013767119015723,  0.10963890178286693708719847752478]
[ -0.82229160187950921144725573239984,  1.0229361963184929117833141996402,  -1.512758814231799805365122632637,   2.3121142197928161050292341106873]
[ -0.48275733697182039898442212778106, 0.62866405653684253473137371430442, 0.70594128160537924726557506770271,  0.14815199882959861698747334577393]
[  0.70334614804537302176239344529766, 0.22388063199966881169984055516601, 0.83071762213161312601662740234805, -0.75794440217665495947886140281076]



</pre>
