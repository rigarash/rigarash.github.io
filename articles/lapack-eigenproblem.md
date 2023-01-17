---
title: "固有値問題におけるLAPACKルーチンの選び方"
emoji: "🧮"
type: "tech"
topics: ["lapack", "eigenproblem"]
published: true
published-at: 2022-12-07
---

## はじめに {#Introduction}

この記事は[数値計算Advent Calendar 2022](https://qiita.com/advent-calendar/2022/numerical_analysis)の7日目の記事です。(Reblog by original author from https://www.rigarash.info/blog/lapack-eigenproblem/ )

```FORTRAN
DSYEV, SSYEVD, CHEEVX, ZHEEVR, SGEEVR...
```

LAPACKには固有値問題を解くための多様なルーチンがあります。過去との互換性をとるためにどんどん複雑になっていますが、実は基本を押さえればどのルーチンを選ぶかは簡単です。

ここでは2022年現在[^1]、固有値問題を解く際のLAPACKルーチン選択の方法を概説します。

**tl;dr ***xyyEVR*** を使えばOK! (x=S,D,C,Z yy=SY,HE,GE)**

基本はこれだけです。ただし、固有値問題の基本とアルゴリズムの背景についても利用者にもわかってほしいので、ざっくりと解説します。

なお、日本には固有値問題の専門家が多数いらっしゃいます[^2]し、日本語文献もたくさんありますので、アルゴリズムや実装に興味がある方はぜひドンドン沼にハマっていってください。

[^1]: 「2022年現在」と書いたもののここ**20年**は変わっていません。

[^2]: [応用数理学会](https://jsiam.org/)に[「行列・固有値問題の解法とその応用」](https://na.cs.tsukuba.ac.jp/mepa/)というその名もズバリな研究部会があります。

## 固有値問題の基本 {#Eigenproblem_basics}

$n \times n$ 行列 $A$ の固有値 $\lambda$ とは、固有方程式 $\| A - \lambda E \| = 0$ の解です。
これは、係数が複素数(または実数)で $\lambda$ の $n$ 次方程式になっていますから、複素数の範囲で重複を含めてちょうど $n$ 個の解をもちます。5次以上の$n$ 次方程式には代数的解法(解の公式)は存在しませんから、これはつまり、固有値は数値的に近似解を求めていく必要があるということです。

LAPACKの実装で注意しなければならないのは、行列の要素が実数か複素数か、という点と、得られた固有値と固有ベクトルの成分が実数か複素数か、が**一致しない**、という点です。表はLAPACKルーチンの1文字目(行列 $A$ の成分)と2,3文字目(行列 $A$ の性質)と固有値/固有ベクトルの関係を示したものです。

| 行列 $A$ の成分 | 行列 $A$ の性質 | 固有値 $\lambda$ の成分 | 固有ベクトル $v$ の成分 |
| ----------- | -------------- | ------ | ------ |
| 実数(S,D)   | 対称(SY)       | 実数   | 実数   |
| 複素数(C,Z) | エルミート(HE) | **実数** | **実数**   |
| 実数(S,D)   | 一般(GE)       | **複素数** | **複素数** |
| 複素数(C,Z) | 一般(GE)       | 複素数 | 複素数 |

行列積のBLASルーチンの選択の際、対称行列に一般行列のアルゴリズムを利用しても計算量が数倍増えるだけですが、固有値問題では実数か複素数かという大きな違いが生じるため、行列の性質を強く意識する必要があります。固有値と固有ベクトルがすべて実数になるのは、実対称行列か複素エルミート行列のときだけです。

### (存在しない) ZSYEV {#ZSYEV}

LAPACKには複素対称行列に対する固有値ルーチンは存在しません。[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)のCMakeLists.txtにはzsyev.fを読み込む部分処理があります[^3]が、これはバグです。複素対称行列は固有値が純虚数になる性質があり、2010年代から理論物理等の分野でたまに出てきますが、LAPACKのAPI設計時にはユースケースがあるとは思われていなかったためか、入っていません。もし複素対称行列の利用が広がってきて、純虚数のデータ構造をどうすべきかといった議論を乗り越えられれば、新規に標準化される可能性もあるのではないか、と思います。

[^3]: EigenがAndroidにも含まれているので検索でたくさんひっかかります。

## 固有値問題のアルゴリズム {#Eigenproblem_algorithm}

ざっくりというと、ハウスホルダー変換、ヘッセンベルグ行列(もしくは実対称か複素エルミートなら三重対角行列)の固有方程式の求解、ハウスホルダー逆変換の流れで行われます。詳しくは(私のような非専門家の解説ではなくて)書籍を読んでください。なお、ハウスホルダー変換などの(ある種)内部ルーチンもLAPACK関数になっているので、固有値ルーチン群は一見数が多いように見えますが、固有値の計算に関しては EV, EVD, EVX, EVR の4種類になります。Driver Routineと呼びます。

### EV/EVD/EVX/EVR の違い {#EV_EVD_EVX_EVR_differences}

三重対角行列の固有方程式の求解部分に差があります。表にまとめると、下のようになります。

| 名称 | アルゴリズム | 演算量 | メモリ量 |
| ---- | ------------ | ------ | -------- |
| EV   | QR           | $6n^3$ | $O(n)$   |
| EVX  | bisection + inverse iteration | $O(n^2) \sim O(n^3)$ | $O(n)$ |
| EVD  | divide & conquer | $O(n^2) \sim O(n^3)$ | $n^2$ |
| EVR  | MRRR         | $O(n^2)$ | $O(n)$ |

アルゴリズム的にも、現実の実装においても、MRRR(えむあーるきゅーびっく)法を用いるのが多くの場合で最速ですし、精度もかなり良いです。MRRR($\textnormal{MR}^3$)法の詳細については、[山本有作先生の解説](https://doi.org/10.11540/jsiamt.15.2_181)を見てください。

実際に、Intel MKLのマニュアルをきちんと読むと、例えばsyevrのreferenceには"[Note that ?syevr is preferable for most cases of real symmetric eigenvalue problems as its underlying algorithm is fast and uses less workspace.](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/lapack-routines/lapack-least-squares-and-eigenvalue-problem/lapack-least-squares-eigenvalue-problem-driver/symmetric-eigenvalue-problems-lapack-driver/syevr.html#syevr)"(EVRを使え)と書いてあります。

## 実際の使い方を学ぶ方法 (#Usage)

これだけでも本の1章分くらいは書けそうですが、ざっくりとまとめると、使い方を学ぶかしかないです。解読しやすいコードとしては、以下のような感じでしょうか。(Fortran95のAPIは使いやすいのでそのまま使えば良いです。)

1. C
   * [Intel MKLのExample](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-lapack-examples/top.html) を参考にする
2. C++
   * 前職で開発したもので、宣伝になりますが、[monolish中のコード](https://github.com/ricosjp/monolish/blob/master/src/internal/lapack/syev/dense_double_syev.cpp) を参考にする
3. Python
   * [scipy.linalg.eigh()](https://github.com/scipy/scipy/blob/v1.9.3/scipy/linalg/_decomp.py#L267-L606) の実装を見てLAPACKの使い方を学ぶ
4. Fortran(95以降)
   * Fortran 95のインターフェースで syevr(a, w) と呼ぶ
