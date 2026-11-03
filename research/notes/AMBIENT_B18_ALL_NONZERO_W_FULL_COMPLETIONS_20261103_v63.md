# 全非零Wへの全完成型の移送

2026-11-03。二つの独立監査済み結果を合わせる、短い系である。

`AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63_REVIEW` は、共通のactual曲線P_Wのalgebraic generic Wで、CA・CB・Nの全完成型を確定した。
`AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63_REVIEW` は、全150係数で

    r^-1 diag(1,1,r,1,1)_* P_W = P_(W/r)

を独立検証した。任意の二つの非零複素値W0,W1に対してr=W0/W1を取ればよい。
generic完成型を持つ非空開集合には非零複素点があるので、これをW0として全W1≠0へ移送できる。
CA・CB・NはそれぞれPGL不変な実被約成分であるから、この移送は個々の成分、その完成環、normalizationとconductorを保つ。

従って**全W≠0**で、CA/CBのprojective完成環は

    C[[g19,s8,t²,s1*t,...,s8*t,t³]] ⊂ C[[g19,s8,t]]

であり、dim28、edim37、depth20、多重度2となる。全81 kernel、conductor、elementary subintegral normalizationはgeneric証明から移る。
Nについては

    C[[g19,s7,t²,s1*t,...,s7*t,t³]] ⊂ C[[g19,s7,t]]

で、dim27、edim35、depth20、多重度2、全64 kernelとなる。
いずれもdomainで、normalizationは正則、finite universal homeomorphismである。

W=0のCA/CB full odd idealやactual CA/CB intersection Dは本系で解決しない。
また、別のV曲線の全特殊値へ移す議論ではない。
再現コードは `code/ambient_B18_all_nonzero_W_completion_20261103_v63.py`。
両reviewのPASSと全入力ハッシュが現在も一致する場合だけ、保存reportを生成する。
