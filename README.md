# 編譯方法
指令：g++ *.cpp -o \<exe\>
# 檔案說明
- maincmd  
  主程式檔案
- algorithm  
  軟判決演算法
- decoding  
  里得所羅門碼的編碼和解碼
- define  
  參數定義檔案
- function  
  計算BER、SER 和 CER 的函數
- Galois  
  加邏瓦域的有限域類別
- Polynomial  
  多項式類別
- global  
  全域(所有檔案)變數的定義
- modulation  
  調變的類別 其中包含兩個類別  
  CModulation 和 CSignal  
  CModulation 為針對 Galois 的調變訊號
  CSignal 則針對 Polynomial
- modulation16PSK  
  繼承 CModulation
- modulation16QAM  
  繼承 CModulation  
- modulation32PSK  
  繼承 CModulation  
- modulation32QAM  
  繼承 CModulation  
- modulationBPSK  
  繼承 CModulation  
  ![image](https://user-images.githubusercontent.com/22883115/136451103-d426ce0a-acbb-4beb-8b2e-6365bad6ff52.png)
