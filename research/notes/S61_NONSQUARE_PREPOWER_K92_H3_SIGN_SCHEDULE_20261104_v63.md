# 実H3係数と外側符号の一意なschedule

2026-11-04。旧guarded source gateは変更しない。そこでのH3_weight_signは
正の記号h lambda/gとの相対符号で、高さ見積もりには影響しない。
四つの未展開H3係数をliteral三行列式として計算する際、そのminusをもう一度掛けてはいけない。

この補助証明はH3六行A,B,D,C,E,Fを形式Laurent monomialで作り、全二十個の三minorを
独立に展開する。十零と十非零の係数・符号を確定し、実scalar program向けの十項scheduleを保存する。
literal四minorには外側Laplace符号だけ、既存六lambdaにはその実係数符号と外側符号の積を使う。
新しいparameter係数や点の計算は行わない。
