# 二H+cubicの残三配置：全六磁気を保持する三つの全一parameter入口

2026-11-08、PRIMARY E、log_manuscript。純係数導出。新operator、minor、point、profile、jet、GBは実行していない。
本稿の新しい二normal formについて別作者review／有限exact readbackは未実施であり、既知親への陽性包含も未証明。
公開130/1084主稿・frontier・全旧producer bytesを変更しない。

## 1. 全源と有限shearの式

`TWO_LINES_PLANE_CUBIC_SUPPORT.md` §1の全源定理を使う。
座標(h,x,y,z,w)、h-Casimir liftで

\[
 P=U\wedge\partial_z+V\wedge\partial_w+(-d w+R)\partial_z\wedge\partial_w,
 \quad [U,V]=dV,
\]

Rはh,x,yの任意二次式である。全六係数を
R=r0 h²+r1 hx+r2 hy+r3 x²+r4 xy+r5 y²とする。
z'=z+ell(h,x,y),w'=w+k(h,x,y)というactual GL5 shearは、baseを保ちdet一、逆はshift符号反転で、

\[
 R\longmapsto R+V(\ell)+(d-U)k.                 \tag{1}
\]

これは有限pushの全式である。base括弧が零なので{ell,k}=0、他の項は{z,k}=−Uk、{ell,w}=Vell及び−d(w'−k)から得る。
infinitesimal像だけを全source正規化と呼んでいない。

## 2. node-flex/tangent-through-flex全源

元既約cubicはC=h y²−x(x−h)²。公開元sourceのinteger scale140を除き、
fiberを(z,w)→(7z/5,20w/7)へactualに変えると

\[
 U=2x(x-h)\partial_x+y(3x-h)\partial_y,\qquad
 V=-2xy\partial_x+[2x(h-x)-3y²]\partial_y,\quad d=x-h.
\]

この時の磁気Rは公開r多項式の四倍である。以下のrはこの正規化後の全六係数を表す。
(1)の全六linear shiftの像は

* ell=h: 0,
* ell=x: −2xy,
* ell=y: 2hx−2x²−3y²,
* k=h: hx−h²,
* k=x: hx−x²,
* k=y: −2xy。

したがって像は正確に四次元で、cokernelはI=r0+r1+r3とJ=r2の二係数。
全域正則な選択

\[
 \ell=(r4/2)x+(r5/3)y,\quad
 k=-(r1+r3)h+(r3-2r5/3)x
\]

によりR'=Ih²+Jhyとなる。rを反転しておらず、全零補集合を保持する。
逆も同じbaseのlinear式の符号反転である。

## 3. cusp-flex/secant-through-flex全源

元cubicはC=h y²−(x+h)³。公開integer scale210を除き、
fiberを(z,w)→(7z/5,30w/7)へactualに変えると

\[
 U=2x(x+h)\partial_x+3xy\partial_y,\qquad
 V=2xy\partial_x+3[y²-(h+x)²]\partial_y,\quad d=x.
\]

磁気Rは公開r多項式の六倍となる。正規化後の全係数について(1)は

* ell=h: 0,
* ell=x: 2xy,
* ell=y: 3y²−3h²−6hx−3x²,
* k=h: hx,
* k=x: −x²−2hx,
* k=y: −2xy。

従って像は四次元、全cokernelはI=r0+r5、J=r2。

\[
 \ell=-(r4/2)x-(r5/3)y,\quad
 k=(2r3-r1)h+(r3+r5)x
\]

で全域R'=Ih²+Jhyとなる。全六磁気の特殊値を除いていない。

## 4. 第一・第二源の全一parameter化と境界

J≠0でfiberを(z'',w'')=(z'/J,w'/J)、tensorをJ倍すると、U,Vと−dwの係数は元に戻り、
磁気だけR→R/Jとなる。よって両源とも全lambda源

\[
 R_\lambda=\lambda h²+hy,\qquad\lambda=I/J\in\mathbb C
\]

へのactual正逆を持つ。公開元のrでのlambdaも順に(r0+r1+r3)/r2、(r0+r5)/r2である。
J=0は元A6の閉補集合として保持する。
有限既知閉和Kへ全lambda源が入ることを別途証明すれば、
PGL/scalar安定性でJ≠0全源がKへ入り、全A6の既約性と閉性からJ=0も全て戻る。
全lambdaの陰性をJ=0へ外挿することは出来ない。

## 5. 第三源との対比と次の正確な有限義務

cusp接線・割線の元全源は既
`TWO_LINES_CUSP_TANGENT_SECANT_MAGNETIC_NORMAL_FORM.md` により
3(Jh²+Iy²)、I=r5、J=r0+r2へ全域actual shearで戻る。
I≠0の全lambda源は3(lambda h²+y²)。その全lambdaがglobalT39/T40外という既算術は受理済みだが、I=0は除外しない。

従って三つの全A6 sourceは、各一つのwhole affine lineとその元閉補集合に正則還元される。
上記二つの新normal formは既固定点の再選択ではなく、全源を保持する有限PGL式である。
既保存D101/M47/source48とも整合し、全源23次元・一般actualorbit22で一parameter商が残ることを説明する。
ただし本稿だけで全parameterのorbit dimension一定又は三sourceがどれかのMPS親へ入るとは結論しない。

次の有限検証が必要なら、元immutable full6 sourceから上記constant frame・全shiftを一回pushし、
全150係数と全逆／全Pfの普遍恒等式だけを確認すればよい。
新M1、Jacobi線形operator、minor、点列、高次normal jetは不要。
そのreadback契約を先にrootへ提出し、既known-parentの指定全曲線が得られるまでprofile探索へ進まない。
