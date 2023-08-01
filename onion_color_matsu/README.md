# onion_color
Onion color estimation by ONT amplicon seq
![概要](https://user-images.githubusercontent.com/108558000/223001174-842e890c-3a36-472a-ba46-8e047b20233f.jpg)


# Primer

|primer name|sequence|
|-----------|--------|
|ANS-F|TTTGCTCGATCGTTTAGCRGAAGAAGA|
|ANS-R|GATCACCATTACACTGATGATGGATC|
|DFR-F|GTTCAACTTTCAAGGGAACTTAGAAATGCC|
|DFR-R_Full_del|TGGGTAGCTATTGGTTCATTCTCTTCA|
|DFR-R_TRN(C→T)|GAGTCGCAACAACGTTAAACGGGTCGT|

※5'末端にバーコード配列を付与したものを用意しておく<br>
バーコード配列は**ONTアンプリコンseq_バーコード配列_1-16.xlsx**の**1stPCR用index配列**シートに記載の<br>
発注配列1と発注配列2を参照してください

# Protocol

## 1st PCR
<img src="https://user-images.githubusercontent.com/108558000/223040836-400b4f6d-5408-4e07-bd6a-2d52d3b146e0.jpg" width="700">


## 2nd PCR
<img src="https://user-images.githubusercontent.com/108558000/223040854-d1ba1d6a-51db-45ac-bf48-8190b013196b.jpg" width="700">


## ゲル抽出
1%アガロースで電気泳動後、1.8~2.1kbp付近のバンドを切り出してゲル抽出をする

## library調製
Ligation Sequencing Kit付属のプロトコルに従ってシーケンス用libraryを調整する<br>
フローセルによって使用するKitのバージョンが異なる
- R9 pore : **SQK-LSK109**または**SQK-LSK110**
- R10 pore : **SQK-LSK112**または**SQK-LSK114**

## MinION flowcellへのアプライ

FLO-MIN114に対して100 fmolをアプライ

# Install
```
$ git clone https://github.com/oku24s/onion_color

$ cp -r /home/usr_share/okumura/onion_color/base_images ./onion_color/container

$ cd onion_color

$ sh container_build.sh
```

※コンテナのベースイメージをlxbioanalysis001の/home/usr_share/okumura/onion_color/base_imagesに保存してある。

container_build.shの内容を変更することでdocker-hub経由で必要なイメージを取得することも可能だが、
container_build.shを実行する前に必要なイメージをコピーしておくことでlxbioanalysis001のローカル環境で全てを完結できる。

# Usage

①configファイルを編集

config.txtを参考に、必要な情報を書き加える

```
indir=/home/okumura/run1/pass                #ONTから出力されたfastqが格納されたディレクトリ

outdir=/home/okumura/run1/output             #出力先のディレクトリ

thread=10                                    #並列処理に使用するスレッド数

sample_list=/home/okumura/run1/sample.csv    #サンプル名を記載したcsvファイル

debarcoding=old_complete                     #debarcodingのための方法 {old, old_complete, v4}

medaka_model=r941_min_hac_g507               #medakaのモデル　{r941_min_hac_g507, r1041_e82_260bps_hac_g632, r1041_e82_400bps_hac_g632}

deepvariant_model=ont_r9_guppy5_sup          #deepvariantのモデル  {ont_r9_guppy5_sup, ont_r10_q20, hifi}

DFR_Full_threshold=10                        #DFR_Fullの足切りリード数

DFR_TRN_threshold=10                         #DFR_TRNの足切りリード数

ANS_h1_threshold=10                          #ANS_h1の足切りリード数

ANS_L_threshold=10                           #ANS_Lの足切りリード数

```

③ パイプラインの実行
```
$ sh ONT_ampliconseq_onion_color.sh
```

↓

対話形式で解析対象とcongig fileの指定を求められるので、configファイルのパスを入力

# その他
DFRがPCRで増幅しない場合、DFR_DTP型またはDFR_LTR型のアリル（いずれも非機能型）の可能性が考えられる。

それぞれ下記のプライマーで検討可能である。
|primer name|sequence|
|-----------|--------|
|DFR-R_DTP|TGAAGCTTGATTAGAGTGCGTTGGAAGATC|
|DFR-R_LTR|CACGTTCCTACTTACTACCAGCGTGCTGAC|
