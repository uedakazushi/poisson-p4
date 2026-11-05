# 同じ固定24の二単元：実sourceと22対22型付き残余

行11を行11−8行44へunimodular変換し、pivot行順[44,11]・列順[49,51]を使う。上三角二対二行列の対角はa=−s0*t0/324000、d=8H*q0/675。実sourceのb=row44,col51もmonomialなので逆は[[1/a,−b/(ad)],[0,1/d]]で、両側全八恒等式を確認する。

残り二十二行・二十二列の各entryは、既存Aij、共有行A44j、既存delta_j=A11j−8A44jの三callbackと小さい単項式係数だけで厳密に指定できる。新しい多項式係数を更新しない。全actual entry・変更しない二十列・旧177行/66制約RHS整合条件を保ち、点ごとのchartは旧底単元・P290・小constraint det・H*s0*q0*t0。固定24のrank24・basisや全逆は未証明。

各future cellについて、既存支持の平行移動の和による支持上界、三つ以下のscalar content分母LCM、primitive整数係数byteからの高さ上界、owned/非圧縮binary上界を記録する。これはSOURCEだけの予算であり、実装・保存addendumや新係数workerのGOではない。多項式GCDや別の行・重み・point・primeは使わない。更なるmonomial構造は一つの既存monomial callbackだけが残る場合に限定して記録し、cancellationをzero-fillしない。
