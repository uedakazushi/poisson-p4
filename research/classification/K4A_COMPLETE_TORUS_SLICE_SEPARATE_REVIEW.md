# K4a complete equivariant slice：別作者 Fraction 全式 review

2026-11-08。**独立 arithmetic／幾何 PASS。**
既存 D101／G48 の全 ambient 証明書を再利用し、同じ canonical point について
新しい polynomial graph の全元式を別実装で確認した。
producer code、SymPy、CAS 又は arithmetic helper を import しない。

## 1. 何を再確認したか

独立 checker は standard-library Fraction と六変数 sparse polynomial
`(c,a,b,u,v,t)` の加法・乗法・微分から Schouten/Euler を再構成した。
新しい点／operator／rank profile／minor の列又は adaptive jet は無い。
既存 D の**同じ**百一 block の determinant を Fraction elimination で再確認し、
既存 all51 action columns の annihilator と N23 の値8/87を直接確認した。
N23 が D kernel に属することも既存 full D に代入した。
全 D/G rank と四十九次元完全 kernel の証明は、先行別作者 ambient review を入力として再利用する。

normalized torus は `s^-1 diag(s²,s,s,1,1)_*`。
全百五十 input weights[−5,2] と三百五十 output characters を定義から計算し、
全 D nonzero support で character が一致することを確認した。
P18 は weight零、N23 は weight−一。
graph の constant は元 P18、linear coefficient は元 N23、
degree二～五 correction は既存百一 raw coordinate columns 内の該当 weight のみである。

producer が保存した百五十×六 polynomial graph をそのまま再構成し、
**全元三百五十 equation、全 degree零～十 coefficient が零**と独立に確認した。
selected百一 equation もその中に含む。
producer の recursive solution を import して結果を繰り返したのではない。

正式 receipt は
`CLASSIFICATION_K4A_COMPLETE_TORUS_SLICE_SEPARATE_REVIEW_V2_20261108.json`。
全照合は0.033秒／13,516 KiBで完了した。
初回 report は job の task_id metadata が inactive E で gate を通らなかったため受理に使わない。
その exact checker/job/report bytes を precontract archive に保持し、
active E_reduced_nonnormal の gate PASS 後、異なる V2 report identity で再実行した。
旧 report を黙って上書きしていない。

さらに actual-curve note の t=174e の整数表を独立に全展開し、
二十五 bracket assignments／三十一 nonzero monomials を saved graph の全百五十係数×五層と比較した。
省略 entry の零も含め全て一致し、全係数が整数である。
この転記確認は新しい Jacobi／rank／profile 計算を含まない。

## 2. 有限 graph が全 formal slice であること

`K4A_EQUIVARIANT_FINITE_SLICE_PROOF.md` の論証を独立に適用する。
既存百一 derivative block が可逆なので selected equations の formal IFT graph は唯一。
slice／selected rows は torus-stable、normal t の weight−一に対して
input weights[−5,2]だから、graph の全 coefficient は degree二～五の一 monomial に限られる。
constant／linear implicit terms は P固定／DN零によって零。
従って今回の literal degree五 family は formal graph の truncation ではなく全 graph。
quadratic full equations の degree十までの全零検証で全 formal residual が零と分かる。

## 3. 正確な local/component 結論

実 Euler-translation／PGL5／scale action と百二次元 slice の map は
G48 と slice tangent の直和から étale。
全 Jacobi scheme を保つ作用なので completed local scheme は smooth四十八変数と
full completed slice C[[t]] の product。
closed orbit 又は Luna slice の仮定は不要。
ゆえに canonical point で raw scheme は smooth四十九次元、
Euler二十四と非零 tensor scale 一を引いた projective scheme は smooth二十四次元。
**canonical orbit closure 二十三次元は full component ではなく、唯一の二十四次元 component 内にある。**

torus-equivariant graph の全 t非零点は actual projective PGL5 同値。
generic nonzero point の orbit tangent は、étale product の四十八方向に
torus が与える t方向を加えた四十九次元。
したがって projective orbit は二十四次元で、この actual orbit closure が上の唯一の component。
この結論にも新しい point-rank profile は用いない。

どの既知成分であるか、prior catalogue と distinct か、成分総数を増やすかは root の別の E義務。
本 review は新 component count を主張しない。

## 4. 再現入口

`classification_k4a_complete_slice_review.py --v2` が正式 checker。
active finite review contract を gate してから実行する。
入力は immutable complete graph／先行 ambient independent certificate／normal readback と
固定した幾何 proof。全実行 checker bytes と hashes は V2 receipt に埋め込む。
現在の producer source を古い証明書の hash へ合わせて書き換える工程は無い。
