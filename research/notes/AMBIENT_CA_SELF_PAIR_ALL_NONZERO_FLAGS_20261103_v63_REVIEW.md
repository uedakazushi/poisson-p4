# 実self-pairの全CA/N fiber：全非零parameterの独立監査

2026-11-03。先行generic finite-open監査は歴史的なscopeのまま保存し、新しい全旗被覆で全a≠0へ強化する。

独立監査済みの75×5 Casimir rank5、175×5 projective factor rank4、三つの全Darboux一次式の分類を再使用する。
全CA proper模型のrank2 preimageはg=0またはf=0で尽くされる。
f=0の全族は二つの定covector Casimirを持つので除外され、残る全sourceはL=<∂y>を持つ。

raw CAの全16 Hamiltonian係数を残した物理式から

    D(ell)=f_k ell/5,  D(k)=−f+f_k k/5

を独立に再検証した。
ell≠xの二candidateはWがxを含むのでk=xと取れる。
5D eigenvalue lambdaに対してf=lambda x²/2+alpha ell²であり、D(x)=−3a x²/5との比較はlambda=2aだけを許す。
従ってell=x+ah−2a²zの既知第二旗だけが残る。
ell=xではD(k)−12a xk/5∈<x²>の全9×4行列を再構成した。
核はspan(h,x)、その二次minorは**50a**なので、単なるQQ(a) generic rankではなく全a≠0で第一旗のWだけを与える。
既知の二つの実150係数sourceが存在するため、全CA geometric fiberはこの二旗で厳密に尽くされる。
固定旗の20係数写像は単射であり、両source微分の単項式minorは全a≠0で非零。
従って全CA scheme fiberは各a≠0で二つの被約点である。

全N sourceは、同じ全flag base上のCA fixed20からC[2]を除いたclosed rank19 subbundleである。
その不変性は先行独立監査で全16 parabolic作用・304列について確認済み。
第一旗のC[2]=0、第二旗のC[2]=−a²なので、全N fiberには第一旗だけが残る。
今回もraw N19と全flag9を物理的に微分し、全150×28と150×9 quotientを再構成した。
9次minorは正確に71663616 a⁹/78125。
従って全N fiberは被約length1であり、proper finite closed-imageとNakayamaから、Nは**全a≠0でcone smooth28、射影smooth27**である。

全CA fiberが各点で二被約点と分かったので、先行のcompletionへの単射・二つのsmooth29 sourceへのsurjectionの証明が全a≠0で適用できる。
reduced CA completionは厳密に二つのsmooth29 formal branchの和集合。
両接空間の和のrank37のminor4076863488 a²⁵/3125は全a≠0で非零なので、actual meeting schemeのedimは21。
0→A→B₁⊕B₂→C→0とdepth lemmaによりcone depth A≤22<29、射影depth≤21<28。
よって個別reduced CAは**全a≠0でnonnormal・singular・非CM**である。

whole CAのembedding dimensionを37と同一視せず、clean meeting、conductor、a=0のfull ringを結論しない。
別の独立監査による同じ点でのCB smoothnessと矛盾しない。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CA_self_pair_all_flags_review_20261103_v63.py --output reports/AMBIENT_CA_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW.json --replay work/v62-review/CA-self-pair-all-nonzero-flags-replay

fresh certificate byte一致、現在hash、全物理行列と厳密minorを検証する。新しいGBは不要。
