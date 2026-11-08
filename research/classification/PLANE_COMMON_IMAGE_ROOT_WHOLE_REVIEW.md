# rank零plane・common-image全四源の独立読解

2026-11-08、root統合担当による純証明の独立確認。
`PLANE_RANK_ZERO_COMMON_IMAGE_WHOLE_RANK_TWO.md` を受理する。
全新operator、rank、minor、kernel又は点の再計算は行わない。

全multiplicity入口・可換common-image二plane・四正常形の全逆は、旧
`PLANE_GLOBAL_COMMONIMAGE_SEPARATE_REVIEW.md` の独立全源定理を用いた。
dependent T,V、scalar、zero symbolをここに含めない。

新しい矢印は、全K1を有限actual shearで除くことである。列の順は
(a_z,b_z,c_z,a_w,b_w,c_w) であり、六shearの逆は符号反転である。
baseにnormal一次式を加える有限代入はweight一を一次にしか変えず、
残りはweight二へ入る。全九P2を同時に動かすため、forcingを捨てる操作ではない。

旧固定matrixを標準ライブラリFractionの直接積だけで読戻した。
SSの全30行は原稿の六shear係数S(s)と一致する。
nilの全30行はS(s)+rNと一致し、差は原kernel第四列だけで、
Nの非零座標はrow10=-1,row14=-1,row21=1である。
旧350行をこのNへ一度掛けても全零である。
この有限恒等式のreceiptは
`reports/CLASSIFICATION_PLANE_COMMON_IMAGE_ROOT_FIXED_MATRIX_READBACK_20261108.json`。
旧D1/D2の新rank又はminorを求めたという意味ではない。

SSでは逆shear後にP1=0なので、旧D2単射性から実P2=0。
nilではN=∂a wedge(-bw∂b-cw∂c+zw∂z)はa非依存のdecomposable tensorで、
strict[N,N]=0。旧D1N=0とD2単射性から、任意のrについて実P2=0となる。
P2の一意性は全actual逆の後に使い、初めからhomogeneous切片へ狭めない。

nilの残tensorにEul wedge(-rw∂a)を加えると、ab/ac成分が消え、
∂a,∂z,∂wだけに支持されることを各項で確認した。
SSも同じ三方向に支持され、いずれもprojective rank二以下である。
raw tensorの一entryからrestriction rankを推定しない。
nilのactual signed Pfは -b(-rcwz²)+c(-rbwz²)=0とも直接確認出来る。

K零二型は旧別作者の全十八higher/全Pf零をそのまま用いる。
k単位の実正逆とk零補集合を合わせ、全四sourceにrank四点は無い。
これは当該common-image laneの全点排除であり、全rank零plane又は全分類の完了ではない。
新成分増分は零。旧producerのlinear-only/all_passed=false等のbytesは変更しない。
