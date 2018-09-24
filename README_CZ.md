# MeshSplit
MeshSplit rozděluje vstupní konvexní model na několik částí/střepů, které dohromady tvoří původní celek.
Jednotlivé střepy jsou tvořeny průnikem původního modelu s buňkami 3D Voroného diagramu sestrojeném nad náhodnými body.

**Vstup:**
- 3d model ve formátu .obj, stačí přetáhnout do okna aplikace

**Výstup:**
- Při uložení výstupu pomocí klávesy "o" ve složce ./out/ jednotlivé kusy 1.obj ... n.obj

**Rozdělení krychle:**

![Cube](https://thumbs.gfycat.com/BrightGentleAmericancurl-size_restricted.gif)

**Rozdělení válce:**

![Cylinder](https://thumbs.gfycat.com/UntimelyThankfulAfricanfisheagle-size_restricted.gif)

**Rozdělení uživatelem zadaného modelu dvacetistěnu:**
![Custom model](https://thumbs.gfycat.com/ImpressionableCooperativeFlounder-size_restricted.gif)

### Jak to funguje

1. Nejprve se sestrojen 3D Voronoi diagram nad náhodnými body. Buňky tohoto diagramu tvorí konvexní obal budoucích střepů.

1. Pro každou buňku:
  Jsou zpracovány trojúhelníky původního modelu. Trojúhelníky ležící mimo buňku jsou zahozeny, trojúhelníky protínající stěnu buňky jsou rozděleny.
  
    - protože stěny buňky dohromady tvoří konvexní tvar, lze test leží/neleží v buňce zjednodušit na test je před/za rovinou tvořící stěnu.
 
1. Zpracované trojúhelníky již tvoří tvar střepu, avšak mají v sobě pořád díry (stěny protínající původní model) které se musí zaplnit. Toto probíhá ve dvou krocích:
 
   1. Výplň stěn obsahující průsečíky:
 
      **Spojení segmentů**
      
        - Segmenty čár vzniklé z průsečíků původních trojúhelníků se stěnou jsou spojeny pokud je to možné.
        
      **Spojení v trojúhelníky**
      
        - Pokud mají dva segmenty stejný krajní bod a stejný vnitřní směr, jsou nahrazeny novým segmentem a trojúhelníkem
        - Zbývající segmenty (pro konvexní mesh) protínají celou stěnu od okraje k okraji
        
       **Spojení posledních segmentů se segmenty stěny**
         - Uzavře poslední místa mezi půnikovými segmenty a segmenty stěny
   1. Výplň stěn bez průsečíků:
      - Stěny, které se celé nacházejí uvnitř modelu nemají žádné průsečíky, ale i tak se musí vyplnit.
      Pokud některý z vertexů prázdné stěny byl vyplněn v předchozím kroku, bude vyplněna i jeho stěna.
      
    ### Omezení
    Program funguje celkem spolehlivě pro konvexní modely. U nekonvexních modelů by bylo třeba vylepšit fázi 3.i) kde vznikají polygony s dírami.
    
    ### Program
    #### Ovládání
    
    Program neobsahuje GUI, ale reaguje na uživatelem stisknuté klávesy:
    
    - **Num 1,2,3** - Přepínání mezi modelem krychle, cylinderu, koule
    - **R** - Restart pro nově vygenerovaný voronoi diagram
    - **W** - Přepínání mezi běžným zobrazením a síťovým modelem (wireframe)
    - **O** - Uložení střepů ve formátu .obj do složky ./out/ v adresáři programu
    
    
    [Download](https://github.com/Fro-Z/MeshSplit/files/2410353/MeshSplit.zip)
      
