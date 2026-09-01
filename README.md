# M1向け RNA-seq配列解析入門

## FASTQから遺伝子ごとのカウントデータを作る

## 1. この演習の目的

RNA-seqでは、シーケンサーから得られた塩基配列を、そのまま遺伝子発現量として利用することはできません。各配列がゲノムのどこに由来するかを推定し、その後、各遺伝子に対応する配列を数える必要があります。

この演習では、paired-end RNA-seqデータから遺伝子ごとのraw countを作成します。

| 段階 | 実行する処理 | 使用するツール | 主な出力 |
| ---: | --- | --- | --- |
| 1 | 配列データの品質を確認する | FastQC | QCレポート（HTML） |
| 2 | アダプター配列と低品質末端を除く | fastp | 前処理済みFASTQ |
| 3 | リードをマウスゲノムへマッピングする | STAR | 座標順ソート済みBAM |
| 4 | BAMファイルを確認する | samtools | BAMの内容、BAMインデックス |
| 5 | リードペアを遺伝子単位で集計する | featureCounts | 遺伝子ごとのraw count |

この演習の到達目標は次のとおりです。

- FASTQ、FASTA、GTF、SAM、BAM、countデータの役割を説明できる
- paired-end FASTQをSTARでマウスゲノムへマッピングできる
- samtoolsを使ってBAMの内容と正常性を確認できる
- featureCountsを使って遺伝子ごとのraw countを作成できる
- 各コマンドの入力ファイルと出力ファイルを区別できる

この演習で扱うのは、**raw countを作成するところまで**です。raw countの読み込み、CPM正規化、低発現遺伝子の除外、発現変化の計算および可視化については、「Pythonデータ解析入門」の`expression_analysis.ipynb`を参照してください。

---

## 2. WSLと作業ディレクトリ

この演習はWSL上で実行します。Windows側のエクスプローラーからWSL内のファイルを開くこともできますが、解析コマンドはWSLのターミナルで実行してください。

### 2.1 作業ディレクトリを作る

```bash
mkdir -p ~/rnaseq_intro
cd ~/rnaseq_intro
```

### 2.2 ディレクトリ構成

```text
rnaseq_intro/
├── data/
│   ├── sample_R1.fastq.gz
│   └── sample_R2.fastq.gz
├── reference/
│   ├── Mus_musculus.GRCm39.dna.primary_assembly.fa
│   ├── Mus_musculus.GRCm39.116.gtf
│   └── star_index_GRCm39_ensembl116/
└── results/
    ├── fastqc/
    │   ├── before_fastp/
    │   └── after_fastp/
    ├── fastp/
    ├── star/
    └── counts/
```

必要なディレクトリを作ります。エクスプローラー上で作成しても構いません。

```bash
mkdir -p data
mkdir -p reference
mkdir -p results/fastqc/before_fastp
mkdir -p results/fastqc/after_fastp
mkdir -p results/fastp
mkdir -p results/star
mkdir -p results/counts
```

教材用FASTQを`data`ディレクトリに配置します。

---

## 3. conda環境とツールのインストール

解析ごとにconda環境を分けると、使用するソフトウェアとその依存関係を整理できます。

### 3.1 conda環境を作る

```bash
conda create \
    -n rnaseq_intro_env \
    -c conda-forge \
    -c bioconda \
    fastqc \
    fastp \
    star \
    samtools \
    subread
```

途中で確認されたら、`y`を入力してEnterを押します。`featureCounts`は`subread`パッケージに含まれています。

### 3.2 conda環境を有効化する

```bash
conda activate rnaseq_intro_env
```

ターミナルの行頭に`(rnaseq_intro_env)`と表示されれば、環境が有効になっています。

### 3.3 インストールを確認する

```bash
fastqc --version
fastp --version
STAR --version
samtools --version
featureCounts -v
```

### 3.4 演習を再開するとき

```bash
cd ~/rnaseq_intro
conda activate rnaseq_intro_env
```

---

## 4. Ensemblからリファレンスをダウンロードする

マッピングにはリファレンスゲノム、遺伝子ごとの集計には遺伝子アノテーションが必要です。今回はマウスの次の組み合わせを使用します。

```text
Assembly: GRCm39
Annotation: Ensembl release 116
```

### 4.1 AssemblyとAnnotation

**Genome assembly**は、染色体の塩基配列と座標を定めたゲノムの版です。同じマウスゲノムでも、assemblyが変わると塩基配列や遺伝子の座標が変わることがあります。今回は`GRCm39`を使用します。

**Gene annotation**は、そのゲノム上のどこに遺伝子、転写産物、exonなどがあるかを記録した情報です。今回はEnsemblがGRCm39に対して作成したrelease 116のannotationを使用します。

例えば、地図にたとえると、assemblyが「地図そのもの」、annotationが「地図上に書き込まれた建物や道路の情報」に相当します。別の版の地図に対するannotationを重ねると座標が対応しないため、マッピングに使うFASTAとカウントに使うGTFは、同じassemblyに対応したものを使用します。

### 4.2 FASTAとは

FASTAは、リファレンスゲノムなどの塩基配列を保存するテキスト形式です。各配列は`>`で始まるヘッダーと、それに続く塩基配列から構成されます。

```text
>1
NNNNNNACGTACGT...
>2
NNNNNNACGTACGT...
```

FASTQと異なり、FASTAには塩基ごとのquality scoreは含まれません。

### 4.3 全ゲノムFASTAをダウンロードする

```bash
cd ~/rnaseq_intro/reference

wget \
https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
```

圧縮ファイルを展開します。

```bash
gzip -d Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
```

FASTAの先頭とヘッダーを確認します。

```bash
head Mus_musculus.GRCm39.dna.primary_assembly.fa
grep '^>' Mus_musculus.GRCm39.dna.primary_assembly.fa
```

識別子は`>`で始まる行の最初の部分のみで、このファイルだと染色体名 (例えば`1`や`2`) が表示されます。後ろの部分には追加情報が含まれていますが、こちらは無視して構いません。

> [!IMPORTANT]
> **問：1番染色体の最初の配列はどうなっていますか？**

<details>
<summary>答を見る</summary>
答：すべて`N`となっています。マウスゲノムの最初の部分がまだ正確に決定されていないためで、`N`は不明な塩基を表します。
</details>

> [!IMPORTANT]
> **問：どのような染色体が含まれていますか？**

<details>
<summary>答を見る</summary>
答：染色体 `1`、`2`、`3`、…など。`MT`, `X`, `Y`はそれぞれミトコンドリア、X染色体、Y染色体を表します。`GL456233`などは未配置の染色体やコンティグを表します。
</details>

### 4.4 GTFとは

GTFは、遺伝子、転写産物、exonなどのゲノム上の位置を記述するタブ区切り形式です。

１列目には染色体、３列目にはアノテーションの種類 (例えば`gene`や`exon`)、４列目には開始位置、５列目には終了位置、７列目にはstrand、が記録されています。9列目には`gene_id`や`transcript_id`などの属性が含まれます。

今回のfeatureCountsでは、主に次の情報を使用します。

- どの領域が`exon`であるか
- 各exonがどの`gene_id`に属するか

### 4.5 全ゲノムGTFをダウンロードする

```bash
wget \
https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz
```

まず、圧縮ファイルを展開します。

```bash
gzip -d Mus_musculus.GRCm39.116.gtf.gz
```

つづいて、先頭を確認します。

```bash
head Mus_musculus.GRCm39.116.gtf
```

先頭にはコメント行が含まれることがあります。

FASTAとGTFで染色体名が一致している必要があります。Ensemblのファイルでは、通常`chr15`ではなく`15`と表記されます。

> [!IMPORTANT]
> **問：GTFファイルの最初の遺伝子のIDや名前は？**

<details>
<summary>答を見る</summary>
答： IDは`ENSMUSG00000104478`、gene_nameは`Gm38212`

9列目に`gene_id "ENSMUSG00000104478"; gene_name "Gm38212"`が含まれています。
</details>

### 4.6 1番染色体用のFASTAとGTFを作る

今回の演習では、計算量を抑えるため、マウス1番染色体だけを含むリファレンスを使用します。

まず、全ゲノムFASTAの検索用インデックスを作ります。

```bash
samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa
```

続いて、配列名が`1`の染色体を取り出します。

```bash
samtools faidx \
    Mus_musculus.GRCm39.dna.primary_assembly.fa \
    1 \
    > Mus_musculus.GRCm39.dna.chromosome.1.fa
```

作成したFASTAのヘッダーを確認します。

```bash
grep '^>' Mus_musculus.GRCm39.dna.chromosome.1.fa
```

出力される配列名が`1`だけであることを確認します。

次に、全ゲノムGTFから、コメント行と1番染色体のannotationを取り出します。

```bash
awk '$1 == "1"' \
    Mus_musculus.GRCm39.116.gtf \
    > Mus_musculus.GRCm39.116.chromosome.1.gtf
```

作成したGTFの先頭を確認します。

```bash
head Mus_musculus.GRCm39.116.chromosome.1.gtf
```

```bash
cd ~/rnaseq_intro
```

---

## 5. FASTQを確認する

### 5.1 今回使用するFASTQ

教員から配布された、次のpaired-end FASTQを使用します。

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

この2ファイルを作業ディレクトリ内の`data`ディレクトリに配置してください。

今回配布するFASTQは、元のRNA-seq FASTQから**リードペアを無作為に10%抽出した教材用データ**です。

ちなみにこのFASTQの元となったファイルは、https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016282/DRX449393/　から入手可能です。

### 5.2 FASTQとは

FASTQは、シーケンサーから得られたリードの塩基配列と、その塩基ごとのquality scoreを保存する形式です。

1リードは4行で表されます。

```text
@read_name
ACGTACGTACGT...
+
FFFFFFFFFFFF...
```

| 行 | 内容 |
| --- | --- |
| 1行目 | リード名 |
| 2行目 | 塩基配列 |
| 3行目 | 区切り記号`+` |
| 4行目 | 各塩基のquality score |

今回のFASTQはgzip形式で圧縮されているため、ファイル名の末尾が`.fastq.gz`になっています。

### 5.3 paired-endとは

RNA-seqでは、cDNAを一定の長さのfragmentにしてから配列を読み取ります。paired-endシーケンスでは、同じfragmentの両端を読みます。

```text
fragment
5' ── R1 → ─────────────── ← R2 ── 3'
```

- `R1`：fragmentの一方の端から読んだ配列
- `R2`：同じfragmentの反対側から読んだ配列

両端の情報を使うことで、片側だけを読むsingle-endよりも、マッピング位置やsplice junctionを推定しやすくなる場合があります。

今回使用するのは1サンプルですが、paired-endなのでFASTQは2ファイルあります。

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

したがって、**1サンプル＝1 FASTQファイルとは限りません**。

### 5.4 ファイルの存在を確認する

```bash
ls -lh data
```

R1とR2の2ファイルがあることを確認します。

### 5.5 gzip圧縮されたFASTQの先頭を表示する

```bash
gzip -dc data/sample_R1.fastq.gz | head
gzip -dc data/sample_R2.fastq.gz | head
```

`head`としているため、それぞれ10行分、すなわち2リード分＋αが表示されます。

### 5.6 リード数を数える

FASTQでは1リードが4行なので、総行数を4で割るとリード数になります。

```bash
gzip -dc data/sample_R1.fastq.gz | awk 'END {print NR / 4}'
gzip -dc data/sample_R2.fastq.gz | awk 'END {print NR / 4}'
```

paired-endでは、R1とR2のリード数が一致することを確認します。異なる場合は、正しいペアとして扱えないリードが含まれている可能性があります。

> [!IMPORTANT]
> 問：R1の最初のリードは何ですか？R1とR2のリード数は？

<details>
<summary>答を見る</summary>
答：CCCAGNTGTCCGAGCTCCTCCTCCCTGTTGAAGCATCTGCACAGGTCCAGCTGTCACTGCTTGGGGACTTGGCCTTのリード、基本的に高品質(E)、一つ低品質(#)な塩基=Nも含まれます。

R1とR2のリード数はいずれも3137561。

</details>

---

## 6. FastQCで配列品質を確認する

マッピング前に、FASTQの品質を確認します。

```bash
fastqc \
    data/sample_R1.fastq.gz \
    data/sample_R2.fastq.gz \
    --outdir results/fastqc/before_fastp
```

```bash
ls -lh results/fastqc/before_fastp
```

R1、R2について、それぞれHTMLファイルとZIPファイルが作成されます。

```text
sample_R1_fastqc.html
sample_R1_fastqc.zip
sample_R2_fastqc.html
sample_R2_fastqc.zip
```

HTMLファイルをWindows側のブラウザで開き、次の4項目を中心に確認します。

### 6.1 Basic Statistics

FASTQ全体の基本情報をまとめた表です。

| 項目 | 確認する内容 |
| --- | --- |
| Filename | 意図したFASTQを解析しているか |
| Total Sequences | リード数。R1とR2で一致しているか |
| Sequence length | リード長。一定長か、長さに幅があるか |
| %GC | GC含量。R1とR2で大きく異なっていないか |

GC含量の適切な値は生物種、発現している遺伝子、ライブラリ調製法によって変わります。特定の値から外れたら直ちに異常というわけではありません。

### 6.2 Per base sequence quality

リードの各位置について、quality scoreの分布を示します。横軸はリード内の塩基位置、縦軸はPhred quality scoreです。

| Phred score | 推定エラー率 |
| ---: | ---: |
| 20 | 1% |
| 30 | 0.1% |
| 40 | 0.01% |

一般に、リードの後半ではqualityが低下しやすい傾向があります。多くの位置でQ30以上であれば高品質と考えられますが、一部の末端でqualityが下がっただけで、直ちにリードを削除する必要があるとは限りません。

R1とR2を比較し、片方だけで著しい低下がないかも確認します。

### 6.3 Per sequence quality scores

各リードの平均quality scoreの分布を示します。

分布のピークが高いquality側にあれば、データの大部分が高品質であることを示します。低quality側に大きなピークや裾がある場合は、低品質なリード集団が含まれている可能性があります。

`Per base sequence quality`が塩基位置ごとの品質を見るのに対し、`Per sequence quality scores`はリード単位の平均品質を見る点が異なります。

### 6.4 Adapter Content

リード内の各位置で、既知のアダプター配列が検出された割合を示します。

アダプター配列は、シーケンサー上でDNA断片を固定したり、PCR増幅したり、サンプルを識別したりするために、ライブラリ調製時に人工的に付加される配列です。生物由来のRNA配列ではありません。

シーケンスするfragmentがリード長より短い場合、fragmentを読み終えた後にアダプター部分まで読み進むことがあります。その場合、リード後半に向かってadapter contentが増加します。

adapter contentがほとんど検出されなければ、通常はアダプター除去の必要性は低いと考えられます。明確に検出される場合は、fastpやCutadaptなどによるアダプター除去を検討します。

### 6.5 `PASS`、`WARN`、`FAIL`の考え方

FastQCの`PASS`、`WARN`、`FAIL`は、一定の基準に基づく機械的な判定です。`WARN`や`FAIL`が1つあるだけで、データ全体が解析不能という意味ではありません。

例えば、RNA-seqでは一部の遺伝子が非常に高く発現するため、配列組成や重複に関する項目で警告が出ることがあります。一方、qualityの大幅な低下や強いアダプター混入は、マッピングに影響する可能性があります。

次の順番で判断します。

1. どの項目に警告が出ているか確認する
2. R1とR2のどちらに出ているか確認する
3. グラフを見て、問題がリード全体か末端だけか確認する
4. ライブラリの種類や解析目的から、マッピングへの影響を考える
5. 必要な場合だけ前処理を行う


> [!IMPORTANT]
> 問：FastQC ReportのPer base sequence qualityとPer base sequence contentを見て、比較的問題のあるリードの位置はどこですか？

<details>
<summary>回答例</summary>
比較的問題のあるリードの位置は、リードの5'末端付近(1-5塩基目)および3'末端（76塩基目）です。
</details>

---

## 7. fastpでFASTQを前処理する

FASTQの品質を確認した上で、必要に応じて前処理を行います。

FASTQの前処理には、アダプター除去や低品質リードのトリミングや除去があります。

アダプター除去は、リードに含まれたアダプター配列を取り除く処理です。低品質トリミングは、主にリード末端のqualityが低い塩基を除く処理です。低品質リードの除去は、全体の平均qualityが低いリードを取り除く処理です。

これらの処理によってマッピングが改善することがありますが、常に行えばよいわけではありません。過剰にトリミングするとリードが短くなり、かえってユニークにマッピングしにくくなる場合があります。

したがって、FastQCで黄色や赤の判定が出たことだけを理由に自動的にトリミングするのではなく、グラフの形と問題の程度を確認して判断します。

### 7.1 fastpの使い方

`fastp`を使って、paired-endの対応関係を保ったまま、アダプター配列の除去と低品質末端のトリミングを行います。

```bash
fastp \
    --in1 data/sample_R1.fastq.gz \
    --in2 data/sample_R2.fastq.gz \
    --out1 results/fastp/sample_R1.trimmed.fastq.gz \
    --out2 results/fastp/sample_R2.trimmed.fastq.gz \
    --detect_adapter_for_pe \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --length_required 30 \
    --html results/fastp/sample_fastp.html \
    --json results/fastp/sample_fastp.json \
    --thread 4
```

| オプション | 意味 |
| --- | --- |
| `--in1`、`--in2` | 処理前のR1とR2 |
| `--out1`、`--out2` | 処理後のR1とR2 |
| `--detect_adapter_for_pe` | paired-endデータについてアダプター配列を自動検出する |
| `--cut_tail` | リード末端から低品質領域を探索して切り取る |
| `--cut_window_size 4` | 4塩基の移動窓でqualityを評価する |
| `--cut_mean_quality 20` | 移動窓の平均qualityが20未満になった末端領域を切り取る |
| `--length_required 30` | 処理後に30塩基未満となったリードペアを除外する |
| `--html`、`--json` | 前処理のレポートを保存する |

`fastp`はpaired-endの2ファイルを同時に処理します。一方のリードが長さなどのフィルターを通過できない場合は、対応する相手側のリードも除かれるため、処理後もR1とR2の対応関係が保たれます。

### 7.2 fastpの結果を確認する

```bash
ls -lh results/fastp
```

HTMLレポートでは、処理前後のリード数、アダプター除去数、quality、リード長などを確認できます。

### 7.3 処理後のFASTQをFastQCで確認する

```bash
fastqc \
    results/fastp/sample_R1.trimmed.fastq.gz \
    results/fastp/sample_R2.trimmed.fastq.gz \
    --outdir results/fastqc/after_fastp
```

処理前後のFastQCを比較し、次の点を確認します。

- `Adapter Content`が減少したか
- リード末端の`Per base sequence quality`が改善したか
- `Total Sequences`がどの程度減少したか
- トリミングによって`Sequence length`に幅が生じたか

前処理後にすべての警告が消える必要はありません。目的は、マッピングを妨げるアダプター配列や低品質末端を減らすことです。

> [!IMPORTANT]
> 問：fastpによって除去されたリードは全体の何パーセントですか？残ったリードの数は？

 <details>
 <summary>回答例</summary>
 約4.4%です。5999054リードが残りました。
</details>

---

## 8. STARインデックスを作成する

続いてSTARによるマッピングを行います。

その前準備として、リファレンスゲノムと遺伝子アノテーションから検索用インデックスを作成します。

今回は、Ensembl release 116のGRCm39の1番染色体のFASTAとGTFから事前に作成した縮小版STARインデックスを作成します。

> [!NOTE]
> 入力FASTQには全染色体由来のリードが含まれていますが、STARが検索するリファレンスは1番染色体だけです。そのため、1番染色体以外に由来するリードはマッピングされません。また、全ゲノムでは複数箇所にマッピングされる配列が、この縮小版リファレンスでは一意にマッピングされたように判定される場合があります。
>
>したがって、今回得られるマッピング率やmulti-mapping率は、通常の全ゲノム解析結果とは直接比較できません。今回は、一連の解析操作を体験するための教材用リファレンスとして使用します。

### 8.1 インデックスの保存先を作る

STARインデックスを保存する空のディレクトリを作ります。

```bash
mkdir -p reference/star_index_GRCm39_chr1_ensembl116
```

### 8.2 STARインデックスを作る

STARインデックスを作成するコマンドは以下の通りです。**1時間以上経っても完了しない場合は、中止して教員に配布インデックスをもらってください。**

```bash
STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir reference/star_index_GRCm39_chr1_ensembl116 \
    --genomeFastaFiles \
        reference/Mus_musculus.GRCm39.dna.chromosome.1.fa \
    --sjdbGTFfile \
        reference/Mus_musculus.GRCm39.116.chromosome.1.gtf \
    --sjdbOverhang 100 \
    --genomeSAindexNbases 12
```

`--sjdbOverhang`は、`リード長 - 1`を指定するのが王道です。今回の場合は、長さが76 bpのリードを想定しているため`75`ですが、デフォルト値の`100`でも実用上問題がないことが多いそうです。

`--genomeSAindexNbases`は、STARの検索用インデックスの内部構造を調整する値です。通常の哺乳類全ゲノムではデフォルトの`14`を使用できますが、1番染色体だけの縮小版リファレンスに合わせて`12`へ小さくします。

ちなみに全ゲノムインデックスは、次のようなコマンドで作成します。20GBを超えるファイルが生成されます。実際の解析では必要ですが、20GBを超えるディスク容量が必要となるため、**こちらは実行しないでください。**

```bash
STAR \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir reference/star_index_GRCm39_ensembl116 \
    --genomeFastaFiles \
        reference/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --sjdbGTFfile \
        reference/Mus_musculus.GRCm39.116.gtf \
    --sjdbOverhang 100
```

### 8.3 作成されたインデックスを確認する

```bash
ls -lh reference/star_index_GRCm39_chr1_ensembl116
```

`Genome`、`SA`、`SAindex`などのファイルが存在することを確認します。

---

## 9. STARでゲノムへマッピングする

STARを使って、R1とR2をマウス1番染色体の縮小版リファレンスへマッピングします。**こちらもかなり時間がかかるかもしれません。3時間以上経っても完了しない場合は、教員に相談してください。**

```bash
STAR \
    --runThreadN 4 \
    --genomeDir reference/star_index_GRCm39_chr1_ensembl116 \
    --readFilesIn \
        results/fastp/sample_R1.trimmed.fastq.gz \
        results/fastp/sample_R2.trimmed.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix results/star/sample_
```

### 9.1 主なオプション

| オプション | 意味 |
| --- | --- |
| `--runThreadN 4` | 4スレッドで実行する |
| `--genomeDir` | STARインデックスの場所 |
| `--readFilesIn` | R1、R2の順にFASTQを指定する |
| `--readFilesCommand zcat` | gzip圧縮FASTQを展開しながら読む |
| `--outSAMtype BAM SortedByCoordinate` | 座標順にソートされたBAMを出力する |
| `--outFileNamePrefix` | 出力ファイル名の先頭部分 |

### 9.2 出力を確認する

```bash
ls -lh results/star
```

| ファイル | 内容 |
| --- | --- |
| `sample_Aligned.sortedByCoord.out.bam` | 座標順にソートされたマッピング結果 |
| `sample_Log.final.out` | マッピング結果の要約 |
| `sample_Log.out` | STAR実行ログ |
| `sample_Log.progress.out` | 実行中の進行状況 |
| `sample_SJ.out.tab` | 検出されたsplice junction |

---

## 10. STARのマッピング結果を確認する

```bash
cat results/star/sample_Log.final.out
```

少なくとも次の項目を確認します。

- Number of input reads
- Uniquely mapped reads number
- Uniquely mapped reads %
- Number of reads mapped to multiple loci
- % of reads unmapped

今回のFASTQは全染色体由来のリードを含みますが、STARインデックスには1番染色体しか含まれていません。そのため、全ゲノムインデックスを使用した場合よりもマッピング率は低くなります。

>[!IMPORTANT]
>問：全部でいくつのリードが入力されたでしょうか。１つの場所のみにマッピングされたリードは何本（何％）でしょうか。複数の場所にマッピングされたリードは何本（何％）でしょうか。
<details>
<summary>答えを見る</summary>
答：全部で2999527×２リードが入力されました。１つの場所のみにマッピングされたリードは482041本（16.07％）です。複数の場所にマッピングされたリードは43302本（1.44％）です。1番染色体のみのインデックスを使用しているため、マッピング率はとても低くなっています。ちなみに、全ゲノムのインデックスを使用した場合、マッピング率は１つの場所のみにマッピングされたリードは2661648本（88.74％）、複数の場所にマッピングされたリードは257828本（8.60％）程度になります。
</details>

---

## 11. samtoolsでSAM/BAMを確認する

### 11.1 SAMとBAM

SAMは、各リードがゲノムのどこにマッピングされたかを記録するテキスト形式です。BAMはSAMとほぼ同じ情報をバイナリ形式で圧縮したものです。

| 形式 | 特徴 |
| --- | --- |
| SAM | テキストなので内容を確認しやすいが、ファイルサイズが大きい |
| BAM | 圧縮されていて小さく、解析ツールで効率よく扱える |

今回はSTARから座標順ソート済みBAMを直接出力しています。SAMを中間ファイルとして保存せず、必要な部分だけ`samtools view`でSAM形式として表示します。

### 11.2 BAMのヘッダーを表示する

```bash
samtools view -H results/star/sample_Aligned.sortedByCoord.out.bam
```

ヘッダーには、リファレンス配列、ソート順、使用したプログラムなどの情報が記録されています。

### 11.3 BAMの先頭をSAM形式で表示する

```bash
samtools view results/star/sample_Aligned.sortedByCoord.out.bam | head
```

| 列 | 名前 | 意味 |
| ---: | --- | --- |
| 1 | QNAME | リード名 (FASTQでのリード名) |
| 2 | FLAG | マッピング状態を表すフラグ|
| 3 | RNAME | マッピング先の染色体 |
| 4 | POS | マッピング開始位置 |
| 5 | MAPQ | mapping quality |
| 6 | CIGAR | リファレンスに対するalignmentの形 |
| 7 | RNEXT | ペアリードの相手の参照名 |
| 8 | PNEXT | ペアリードの相手のマッピング開始位置 |
| 9 | TLEN | テンプレート長 （1つのリードの先頭からもう1つのリードの末尾まで）|
| 10 | SEQ | リード配列 |
| 11 | QUAL | クオリティスコア （FASTQでの品質スコア）|

#### 11.3.1 FLAG

FLAGは、各リードのマッピング状態を表したものです。

| 数字 | 意味 |
| ---: | --- |
| 1 | ペアエンドリード |
| 2 | マッピング成功（ペアリードも） |
| 4 | マッピング失敗 |
| 8 | ペアリードがマッピング失敗 |
| 16 | リバース鎖マッピング |
| 32 | ペアリードがリバース鎖マッピング |
| 64 | リード１(どっちのFASTQ由来か)|
| 128 | リード２ （どっちのFASTQ由来か）|
| 256 | ベストじゃないアラインメント（multiple mappingしたときの） |

マッピング成功したペアリードはそれぞれ「99と147(1+2+32or16+64or128)」あるいは「83と163(1+2+16or32+64or128)」となります。

片方しかマップされなかったペアリードは、それぞれ「73(1+8+64)と133(1+4+128)」などになります。

256より大きい場合は、複数の位置でマッピングするリードのベストではないアラインメントとなります。複数の位置でマッピングするリードのベストのアラインメントは、ユニークな位置にマッピングされたものと同じFLAGになるため、**FLAGだけではマルチプルマッピングかどうかは断定できず、MAPQを見る必要があります**。

#### 11.3.2 MAPQ

MAPQは、mapping qualityはアラインメントの信頼度を表しますが、ソフトウェアによってやや異なった値が用いられるので注意しましょう。STARで覚えておくべきのは以下の3種類のみです。

| MAPQ | 意味 |
| --- | --- |
| 255 | 品質が計算されていない、ユニークなマッピングを含む |
| 1-10 | マッピング品質が低い (信頼度が低い) |
| 0 | 複数箇所に全く同じスコアでマッピング、リピート配列由来が多い |

ちなみに、エピゲノム解析で使うこともあるBowtie2では、ユニークなマッピングに相当するのは42です（なんで？）。

#### 11.3.3 その他の列

CIGARは、リファレンスに対するalignmentの形を表します。例えば、`10M1I5M`は、10塩基がマッチし、1塩基が挿入され、さらに5塩基がマッチしていることを示します。`7S68M`のように、`S`はソフトクリップ＝リードの一部がリファレンスにマッチせずに切り取られていることを表します。

RNA-seqでは、exonをまたぐリードが存在します。CIGARに含まれる`N`は、リファレンス上で読み飛ばされた領域を表し、splice junctionに対応する場合があります。

TLENは、テンプレート長を表します。ペアリードの場合、1つのリードの先頭からもう1つのリードの末尾までの距離です。

QUALは、リード配列のクオリティスコアを表します。

> [!IMPORTANT]
> 問：２番目のアラインメントのCIGARは何でしょうか。何を意味しているでしょうか。

<details>
<summary>答えを見る</summary>
答：２番目のアラインメントのCIGARは`7S68M`です。これは、最初の7塩基がソフトクリップされ、続く68塩基がリファレンスにマッチしていることを意味しています。
</details>

### 11.4 BAMファイルが正常に作成されたことを確認する

STARによって作成されたBAMファイルを確認します。

```bash
samtools quickcheck -v results/star/sample_Aligned.sortedByCoord.out.bam
```

`samtools quickcheck`は、BAMファイルのヘッダーや末尾を調べ、ファイルが途中で切れていないかを簡易的に確認するコマンドです。問題がなければ、通常は何も表示されません。

BAMに含まれる**アラインメント**の内訳は次のコマンドで確認できます。

```bash
samtools flagstat results/star/sample_Aligned.sortedByCoord.out.bam
```

STARのログは主に**リードペア単位**、`flagstat`はBAM内のR1・R2それぞれのrecord単位で集計します。そのため、`flagstat`の`mapped`の割合をFASTQ全体のマッピング率として解釈することはできません。

---

## 12. BAMインデックスを作る

BAMインデックスは、BAM内の特定のゲノム領域へ高速にアクセスするための索引です。

```bash
samtools index results/star/sample_Aligned.sortedByCoord.out.bam
ls -lh results/star/sample_Aligned.sortedByCoord.out.bam* 
```

通常、BAMと`.bam.bai`が表示されます。

### 12.1 特定領域のalignmentを取り出す

例として、各染色体にマッピングされたリード数を確認します。

```bash
samtools idxstats results/star/sample_Aligned.sortedByCoord.out.bam | head
```

続いて、任意の領域を指定してalignmentを表示できます。次は1番染色体の一部を指定する例です。

```bash
samtools view results/star/sample_Aligned.sortedByCoord.out.bam 1:30000000-30100000 | head
```

このような領域指定では、BAMが座標順にソートされ、インデックスが作成されている必要があります。指定した領域にリードがなければ、何も表示されません。

> [!NOTE]
> BAMインデックスはfeatureCountsの実行には必須ではありません。特定領域の抽出やIGVでの表示など、BAMへランダムアクセスするときに必要です。

---

## 13. featureCountsで遺伝子ごとのcountを作る

featureCountsは、BAMに記録されたalignmentと、GTFに記録されたexonの位置を重ね合わせ、各遺伝子に割り当てられるfragment数を集計します。

```bash
featureCounts \
    -a reference/Mus_musculus.GRCm39.116.gtf \
    -o results/counts/gene_counts.tsv \
    -t exon \
    -g gene_id \
    -p \
    --countReadPairs \
    -s 2 \
    results/star/sample_Aligned.sortedByCoord.out.bam
```

### 13.1 主なオプション

| オプション | 意味 |
| --- | --- |
| `-a` | 遺伝子アノテーションとして使用するGTF |
| `-o` | count結果の出力先 |
| `-t exon` | GTFの`exon`を対象とする |
| `-g gene_id` | 同じ`gene_id`を持つexonを遺伝子単位にまとめる |
| `-p` | paired-endデータとして扱う |
| `--countReadPairs` | readではなくfragment単位で数える |
| `-s 2` | reverse-strandedライブラリとして数える |

今回のデータはTruSeq Stranded mRNA Library Prepで作成されたreverse-strandedライブラリなので、`-s 2`を指定します。strandednessが異なるデータでは、この設定も変更する必要があります。

paired-endではR1とR2が同じfragmentに由来します。`-p --countReadPairs`を指定することで、R1とR2を別々の2リードではなく、1 fragmentとして数えます。

---

## 14. count結果を確認する

### 14.1 countデータとは

countデータは、各遺伝子に割り当てられたRNA-seq fragment数を表す表形式データです。featureCountsで得られる値は、正規化前の**raw count**です。

### 14.2 出力ファイルを確認する

```bash
ls -lh results/counts
```

主に次の2ファイルが作成されます。

```text
gene_counts.tsv
gene_counts.tsv.summary
```

### 14.3 countファイルの先頭を表示する

```bash
head results/counts/gene_counts.tsv
```

| 列 | 内容 |
| --- | --- |
| Geneid | Ensembl gene ID |
| Chr | 染色体 |
| Start | 集計対象領域の開始位置 |
| End | 集計対象領域の終了位置 |
| Strand | 遺伝子のstrand |
| Length | 集計に用いられたexon領域の長さ |
| BAMファイル名の列 | その遺伝子に割り当てられたraw count |

> [!IMPORTANT]
> 問：遺伝子PTMA (ENSMUSG00000026238)のraw countはいくつか？ ※grepで確認してみましょう。

<details>
<summary>答えを見る</summary>
答：合計8181フラグメントが遺伝子PTMA (ENSMUSG00000026238)に割り当てられています。
</details>

### 14.4 assignment summaryを確認する

```bash
cat results/counts/gene_counts.tsv.summary
```

summaryでは、それぞれのfragmentが遺伝子へ割り当てられたか、割り当てられなかった場合は何が理由だったかを確認できます。

| 項目 | 意味 |
| --- | --- |
| `Assigned` | 遺伝子に割り当てられたfragment |
| `Unassigned_NoFeatures` | 対象とするexonに重ならなかったfragment |
| `Unassigned_Ambiguity` | 複数の遺伝子に重なり、割り当て先を1つに決められなかったfragment |
| `Unassigned_MultiMapping` | ゲノム上の複数箇所にマッピングされたfragment |
| `Unassigned_Unmapped` | マッピングされていないfragment |

STARの`Uniquely mapped reads number`は、ゲノム上の1か所にマッピングされたリードペア数です。一方、featureCountsの`Assigned`は、そのうちGTFに記載された遺伝子のexonへ割り当てられたfragment数です。

今回のfeatureCountsのオプションでは、次の関係がおおむね成り立ちます。

```text
STARのUniquely mapped reads number
    ≒ Assigned
    + Unassigned_NoFeatures (エクソン以外の領域へのマッピング)
    + Unassigned_Ambiguity (重複する遺伝子の重複部分へのマッピング)
```

一意にマッピングされたfragmentであっても、イントロン、遺伝子間領域、GTFに記載されていない領域にマッピングされた場合は、遺伝子のcountには入りません。また、今回のように`-s 2`を指定した場合、遺伝子のstrandと期待される向きが一致しないfragmentも、主に`Unassigned_NoFeatures`に含まれます。

完全な等式ではなく`≒`としているのは、ソフトウェア間で判定と集計の方法が異なり、追加のフィルターを指定した場合には、ほかの`Unassigned_*`項目も考慮する必要があるためです。

> [!NOTE]
> もしマルチプルマッピングを配分したcountを計算したい場合は、別のカウントソフトウェア（16.4参照）を使用する必要があります。

> [!IMPORTANT]
> 問：遺伝子に割り当てられたリードはいくつか？最初のFASTQファイルのうちの何割になったか？

<details>
<summary>答えを見る</summary>
答：遺伝子に割り当てられたリード数は407380です。元のリード数が3137561なので、最初のFASTQファイルのうちの約13%が遺伝子に割り当てられたことになります。
</details>

---

## 15. 次の解析への接続

本演習では、1サンプルのpaired-end FASTQからraw countを作成しました。

複数サンプルのraw countがそろうと、サンプル間の発現量や発現変化を解析できます。続きは、「Pythonデータ解析入門」の`expression_analysis.ipynb`を参照してください。

そちらのNotebookでは、主に次の処理を扱います。

- countデータの読み込み
- サンプルごとの総カウントの確認
- CPM正規化
- 低発現遺伝子のフィルタリング
- `log2(CPM + 1)`変換
- 発現パターンのヒートマップ
- log2 fold changeの計算
- 発現増加・減少遺伝子の抽出
- 解析結果の保存

今回得られた1サンプルのcountだけでは、群間比較や発現変動解析はできません。発現変動を評価するには、比較対象となる複数の条件と、生物学的反復を含む複数サンプルが必要です。

---

## 16. 今回は実行しなかった工程

### 16.1 MultiQC

MultiQCは、複数サンプルのFastQCやマッピング結果を1つのレポートにまとめるツールです。今回は1サンプルだけなので実行しません。多数のサンプルを扱う場合は、サンプル間で品質やマッピング率を比較するために有用です。

### 16.2 rRNAや汚染配列の評価

mRNA-seqでは、ライブラリ調製時にpoly(A) RNAを選択したり、rRNAを除去したりします。この処理が不十分だと、rRNA由来のリードが多くなり、遺伝子発現解析に利用できるリードの割合が低下します。

rRNA混入は、rRNA配列へのマッピング、FastQCの`Overrepresented sequences`、専用のQCツールなどで評価できます。別生物種やアダプターなどの汚染については、複数の候補ゲノム・配列へ照合する方法があります。

今回は解析の中心をFASTQからgene countまでに絞るため、これらの評価は実行しません。

### 16.3 duplicateの評価・処理

duplicateは、同じ位置へ同じようにマッピングされたリードまたはfragmentです。PCR増幅によって生じるtechnical duplicateもありますが、RNA-seqでは高発現遺伝子から同じfragmentが独立に多数得られることもあります。

そのため、通常のbulk RNA-seqでduplicateを一律に除去すると、本当に高発現している遺伝子のカウントまで減らす可能性があります。一般的な遺伝子発現解析では、duplicateを機械的に除去しません。

UMIを付加したライブラリでは、UMIを利用してPCR duplicateを識別できるため、扱いが異なります。

### 16.4 transcript-level定量

featureCountsは、今回はexonに重なるfragmentを`gene_id`単位で集計しています。一方、1つの遺伝子から複数の転写産物が作られる場合、どのisoformに由来するリードかを区別したいことがあります。

RSEM、Salmon、kallistoなどは、複数の転写産物に対応し得るリードを統計的に割り当て、transcript-levelの推定カウントやTPMなどを計算します。

transcript-level定量では、共有exonに由来するリードをどの転写産物へ割り当てるかという不確実性があります。今回はSAM/BAMとgene-level countの関係を学ぶことを優先し、扱いません。

### 16.5 発現量の正規化と発現変動解析

CPMなどへの正規化、低発現遺伝子の除外、発現変化の計算および可視化は、「Pythonデータ解析入門」を参照してください。

### 16.6 workflow化

実際の解析では、複数サンプルに対して同じコマンドを実行するため、シェルスクリプト、Snakemake、Nextflowなどを使って処理を自動化することがあります。プログラミングのように解析フローをまとめて管理できるので、実際に解析を行う際には１ファイルにフローをまとめることも検討してみましょう。

---

## 17. 参考資料

- [STAR公式リポジトリ](https://github.com/alexdobin/STAR)
- [SAMtools公式ドキュメント](https://www.htslib.org/doc/)
- [featureCounts公式ページ](https://subread.sourceforge.net/featureCounts.html)
- [FastQC公式ページ](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [fastp公式リポジトリ](https://github.com/OpenGene/fastp)
- [Ensembl FTP Download](https://www.ensembl.org/info/data/ftp/index.html)
- [Ensembl Mus musculus](https://www.ensembl.org/Mus_musculus/Info/Index)
- [Bioconda: STAR](https://bioconda.github.io/recipes/star/README.html)
- [Bioconda: samtools](https://bioconda.github.io/recipes/samtools/README.html)
- [Bioconda: subread](https://bioconda.github.io/recipes/subread/README.html)
- [Bioconda: FastQC](https://bioconda.github.io/recipes/fastqc/README.html)
- [Bioconda: fastp](https://bioconda.github.io/recipes/fastp/README.html)
