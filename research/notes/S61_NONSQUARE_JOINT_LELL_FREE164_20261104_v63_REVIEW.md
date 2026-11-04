# actual free164→122：117 pivot restoration の独立監査

2026-11-04。対象は凍結済み `S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63`。
既存 mixed free183 から単純に19行を引く論法ではなく、元319 sourceでquadratic117 pivotを
先に行い、その後既知38方向を具体的に移す有限circuitを監査する。
全164×122 coefficient matrix、rank、新点、素数、GBは計算しない。

## 1. 安全なsource/target順序

受理済みunit98にoverlap5とqpair14を加えた117 source/target indexを独立に再構成する。
全pivot sourceは元prime quadratic174にあり、その補集合57と合わせて174をちょうど分割する。
targetも117と122に分割する。従ってfull319は、同じcokernelの

    quadratic57 + nonquadratic145 → target122

という202→122表示へ移る。非二次145は元cubic112とquartic33のlabelのままである。
全33quarticをこの段階で保持する。

## 2. 係数の逆輸送を検証する

forward generator操作が g'_r=g_r/p、g'_t=g_t−c_t g'_r のとき、
係数の逆変換は

    k_r=(k'_r−Σ_t c_t k'_t)/p,
    k_t=k'_t (t≠r)

である。この転置・符号・pの位置を固定して、保存した14段を逆順、次いで5段を逆順に戻す。
各段のpivotと全elimination coefficientsは直前のindependent QPAIR/OVERLAP REVIEWに束縛される。
overlap5ではqpair58行への係数も必要であるため、full NF98 operatorと元174 quad productsから
**全58×5 scalar roots**を再構成し、保存値と照合する。
genuine cubic/quarticの係数は一つも展開しない。

初期57係数を未選択sourceに置き、19選択source係数を零として逆操作を合成する。
その76 outputを既存unit98 restorationへ入れ、prime174、old174、oldP22 output22の全rootを確認する。
これが実際のP57である。K57の全57rootは旧K76の未選択57座標への制限として確認する。

ここで「未選択係数は変わらない」は単なる推測ではない。
各操作はpivot generatorだけを正規化し、他行からそれを引くので、逆係数式ではpivot係数だけが変わる。
一度もpivotにならない57座標は全19段を通じて不変である。
同じ理由で非二次source145の係数も不変なので、旧15のcubic/quartic unit pivotはliteralに残る。

## 3. P57K57=Iのexact proof

旧K22は元pure quadratic sourceのliteral kernelである。
117 unit pivots後、selected sourceとselected targetのsquareはunit triangular。
したがって変換されたpure kernelの117 selected source係数は全て零となる。
57座標へ制限して上の逆係数操作と98restorationを行えば、旧K174を一意に復元する。
受理済みoldP22*oldK22=I22から、**P57*K57=I22**。

この議論はpure kernel上のexact identityとunit squareだけを使い、generic rankから全点へ移す推論、
quadratic cokernel flatness、quad57のunimodular minor、projective quotientの自由性を仮定しない。
一般のfull kernelの二次部分を単独で復元して元に戻る、という不要な仮定もしない。

## 4. 既知38方向の補空間を構成する

旧15を不変な非二次pivotで消すとsource187=quad57+surplus130。
Bnew22はこの **非二次surplus130** から選び、separatorの元319座標182,186を避ける。
残りC108に両座標が残ることを確認する。quad57からB22を選んでいない。

P57K57=Iによる source可逆分解

    k=P57*x−u,
    z=x−K57*k,
    (x,u,h)=(K57*k+z,P57*z,h)

により自由165表示 [A57+Bnew*P57,C108] を得る。
新kernel165の全rootも、実kernelの不変57係数、旧u22/h108、P57/K57から再構成する。
separatorは座標57,58にあり、その和は既知unit theta=−1024000aC³J/243。
これは既に受理したfree281/free183と同じscalar identity rootそのものである。
定数座標和変換後にこのunit kernel成分を削除して、自由164→自由122の表示が完成する。
従ってこれは元319と同じactual cokernelを持つ。旧sourceの固定164行部分行列だとは主張しない。

## 5. 独立circuit検証とscope

保存scalar22244 node、linear4616 nodeを別のtyped expression internerで読み、
全復元root、P57/K57、実kernel165、target141→136→122およびfull239→122を独立に構成して照合する。
coordinate/old98/P22のprefixは受理済みfree183とbyte-levelで一致する。
nativeは全1140 atomの分母が既知a,d,F,B,C,Jの積であることと、同じ小scalar separator identityだけを検証する。
大きいrootはfactored expressionのままである。

known kernel directions38、残り42は変わらない。
free164表示が完成してもtarget122の全射性、full239、全kernel exhaustion、全点cubic57は未証明。

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_joint_Lell_free164_review_20261104_v63.py --output work/marked-square-v63/joint-free164-review-REPLAY
```
