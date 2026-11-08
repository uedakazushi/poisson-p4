# Ga9 linear-w：有限全支持被覆からの23次元成分性

2026-11-08、PRIMARY E、全既知親非包含と既受理の有限閉被覆を結ぶ純粋な consequence。新しい point、ambient operator、rank/minor、profile 又は jet は無し。catalogue／成分数／主稿／Gitはrootのみが更新する。

## 入力を分離する

1. `NIL211_SHORT_TRANSLATION_REVIEW.md` と `GA9_LINEAR_W_FIXED_DIVISOR_MAGNETIC_NORMAL_FORM.md`：linear-wの全十一sourceのPGL飽和閉包 C_G は、一つのactual23軌道 O の閉包。全磁気零境界も同じ閉包に含む。canonical actual divisor は −4h(4h³w−8h²y²+y⁴)、被約H＋絶対既約quartic。
2. `GA9_LINEAR_W_PRIOR142_WHOLE_COMPONENT_EXCLUSIONS.md` と factor-image scope、最終142-row receipt：O は既知142成分の閉和 K の外。T37の新しいwhole simple-component低次数surface定理とGa9 degree九の別作者受理を含む。
3. `KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md` と root 独立 `KNOWN142_SEVEN_CANDIDATE_GLOBAL_COVER_SEPARATE_REVIEW.md`：全被約Poisson支持 P_red は

```
K ∪ C_G ∪ C_N ∪ C_F ∪ C_C ∪ C_D ∪ C_T ∪ C_U
```

という有限既約閉和。plane三候補は元全六磁気sourceの閉飽和、constantdouble/tripleとunimodularは各全源を含む軌道閉包である。一般点だけのG上界からこの全支持等式を仮定したのではない。

## 他の六候補は閉factor像だけで除ける

全因子collisionを含む proper multiplication image を使う。

* C_N,C_F,C_C の全actual一般式は二つのsimple線形因子とcubicの積。従って全閉包は M113={L1*L2*G3} の affine cone内（Pfaffian零も含む）。
* C_D,C_T はactual二重Hとcubic、従って全閉包は M23^rep={L²*G3} 内。
* C_U はactual5H、従って全閉包は M5={L⁵} 内。

Ga9の非零canonical divisorは単純Hと既約quarticだけなので、この三閉像のいずれにも入らない。M113でL1=L2などのcollisionを許しても、第二のlinear factor又は重複Hを必要とすることは変わらない。irreducibilityやreducednessを閉条件とは呼ばない。

よって任意の p∈O は、全被覆の他の全七項 K,C_N,C_F,C_C,C_D,C_T,C_U の外にあり、唯一 C_G の項だけに属する。

## 成分性の有限閉和論法

pを含むP_redの既約成分 Z を取る。PGL作用は連結なので各既約成分を保存し、ZはOとその閉包C_Gを含む。一方、有限閉和を覆う既約Zはその一つの閉項に含まれる。pは他の全項の外だから Z⊂C_G。従って

```
Z=C_G=closure(O),   dim C_G=23.
```

これで全linear-w源は一つの新しい全空間被約支持成分に入る。別の六候補の相互成分性・相違を先に証明する必要はない。dense軌道の外の磁気零境界が旧成分と交わることも、この論理と両立する。

この議論は full Jacobi scheme が smooth／reduced であるとは主張しない。canonical tangent dimension、nilpotent厚み、normal sliceの完了は必要でなく、計算もしない。exact全成分数は他の有限候補の判定とroot台帳登録に依存する。
