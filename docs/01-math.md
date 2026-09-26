# 1. Математика и пространственные структуры

[← Оглавление](README.md) · [Твёрдые тела →](02-rigid-bodies.md)

**Что это и зачем.** `src/math/` — небольшая библиотека линейной алгебры без шаблонов: `float`-векторы и матрицы для горячих циклов физики, `double`-матрица произвольного размера для точных решений. `src/spatial/` — две иерархии ограничивающих объёмов: статическая BVH (строится один раз, для мешей) и динамическое AABB-дерево (обновляется по мере движения тел).

Всё подключается одним заголовком [src/math/Math.h](../src/math/Math.h).

| Файл | Типы и функции |
|---|---|
| [Scalar.h](../src/math/Scalar.h) | `kPi`, `kInf`, `clampv`, `sqr`, `lerp`, `degToRad` |
| [Vector2.h](../src/math/Vector2.h), [Vector3.h](../src/math/Vector3.h), [Vector4.h](../src/math/Vector4.h) | векторы, `dot`, `cross`, `normalize`, `anyPerpendicular` |
| [Quaternion.h](../src/math/Quaternion.h) | `Quaternion`, `slerp`, `extractRotation` |
| [Matrix3x3.h](../src/math/Matrix3x3.h) | `Matrix3x3`, `symmetricEigen` (Якоби 3×3) |
| [Matrix4x4.h](../src/math/Matrix4x4.h) | аффинные преобразования, `perspective`, `lookAt`, `inverse` |
| [MatrixNxN.h](../src/math/MatrixNxN.h) | `MatrixNxN` (LU, Холецкий, Якоби), `solveSmall` (n ≤ 4) |
| [AABB.h](../src/math/AABB.h) | ограничивающий параллелепипед, slab-тест луча |

---

## 1.1 Векторы

`Vector3` — три `float` без выравнивания и SIMD-магии. Все операции `constexpr` или `inline`. Полезные нестандартные функции:

- `anyPerpendicular(a)` — какой-нибудь единичный вектор, перпендикулярный `a`. Используется для построения касательного базиса контакта. Выбирается по величине `a.x`, чтобы не делить на почти ноль:

[src/math/Vector3.h:49](../src/math/Vector3.h#L49)
```cpp
inline Vector3 anyPerpendicular(const Vector3& a) {
    return normalize(std::fabs(a.x) > 0.57f ? Vector3(a.y, -a.x, 0) : Vector3(0, a.z, -a.y));
}
```

Порог $0.57 \approx 1/\sqrt 3$: хотя бы одна компонента единичного вектора не меньше этого значения, поэтому длина результата всегда отделена от нуля.

`normalize` возвращает нулевой вектор для вырожденного входа (длина < 1e-20) — вызывающий код не получает NaN.

---

## 1.2 Кватернионы

Ориентация хранится единичным кватернионом

$$
q = (w, x, y, z) = \left(\cos\tfrac{\theta}{2},\ \mathbf n \sin\tfrac{\theta}{2}\right),
$$

где $\mathbf n$ — ось поворота, $\theta$ — угол. Поворот вектора: $\mathbf v' = q\,\mathbf v\,q^{-1}$; в коде — через матрицу `toMatrix3x3()`.

### Интегрирование угловой скорости (экспоненциальное отображение)

За шаг $\Delta t$ тело с угловой скоростью $\boldsymbol\omega$ поворачивается точно на угол $|\boldsymbol\omega|\Delta t$:

$$
q(t+\Delta t) = \exp\!\left(\tfrac12 \boldsymbol\omega \Delta t\right) q(t), \qquad
\exp(\mathbf r) = \left(\cos|\mathbf r|,\ \frac{\sin|\mathbf r|}{|\mathbf r|}\,\mathbf r\right).
$$

При малом угле отношение $\sin\theta/\theta$ даёт $0/0$. Код переходит на ряды Тейлора при $\theta < 10^{-2}$:

[src/math/Quaternion.h:96](../src/math/Quaternion.h#L96)
```cpp
Quaternion integrated(const Vector3& omega, float dt) const {
    const float th = 0.5f * dt * std::sqrt(omega.x * omega.x + omega.y * omega.y + omega.z * omega.z);
    float c, sinc; // cos(th), sin(th)/th
    if (th < 1e-2f) {
        const float t2 = th * th;
        c = 1.0f - t2 * (0.5f - t2 / 24.0f);           // 1 - th^2/2 + th^4/24
        sinc = 1.0f - t2 * (1.0f / 6.0f - t2 / 120.0f); // 1 - th^2/6 + th^4/120
    } else {
        c = std::cos(th);
        sinc = std::sin(th) / th;
    }
    const float k = 0.5f * dt * sinc; // vector part = omega/|omega| * sin(th) = omega * dt/2 * sinc
    Quaternion e{c, omega.x * k, omega.y * k, omega.z * k};
    return (e * (*this)).normalized();
}
```

Обратная операция — **логарифм** `log()` ([Quaternion.h:80](../src/math/Quaternion.h#L80)): вектор поворота $\mathbf r = \theta\,\mathbf n$ по кратчайшей дуге (при $w < 0$ кватернион сначала меняет знак). Логарифм нужен везде, где ошибку ориентации надо превратить в вектор: сочленения, блокировка вращения контактов, CCD.

### Сферическая интерполяция (slerp)

$$
\operatorname{slerp}(a, b, t) = \frac{\sin((1-t)\Omega)}{\sin\Omega}\,a + \frac{\sin(t\Omega)}{\sin\Omega}\,b, \qquad \cos\Omega = a\cdot b.
$$

[src/math/Quaternion.h:116](../src/math/Quaternion.h#L116)
```cpp
inline Quaternion slerp(const Quaternion& a, Quaternion b, float t) {
    float c = dot(a, b);
    if (c < 0) { b = Quaternion{-b.w, -b.x, -b.y, -b.z}; c = -c; } // same rotation, shorter way round
    if (c > 0.9995f) { // nearly parallel: linear interpolation is exact enough and stable
        Quaternion r{a.w + (b.w - a.w) * t, a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t};
        return r.normalized();
    }
    float theta = std::acos(c);
    float wa = std::sin((1 - t) * theta) / std::sin(theta), wb = std::sin(t * theta) / std::sin(theta);
    return Quaternion{a.w * wa + b.w * wb, a.x * wa + b.x * wb, a.y * wa + b.y * wb, a.z * wa + b.z * wb}.normalized();
}
```

Две детали: $q$ и $-q$ — один и тот же поворот, поэтому при $a\cdot b < 0$ берётся $-b$ (короткий путь); при почти совпадающих кватернионах $\sin\Omega \to 0$ и используется нормированная линейная интерполяция.

### Вращательная часть матрицы: `extractRotation` (Müller et al. 2016)

Задача: дана матрица деформации $\mathbf A$ (например, $\sum_i (\mathbf p_i - \mathbf c)\,\mathbf q_i^{\mathsf T}$ в shape matching, гл. 3). Нужна ближайшая к ней матрица поворота $\mathbf R$ из полярного разложения $\mathbf A = \mathbf R\mathbf S$.

Классический путь — SVD или $\mathbf R = \mathbf A(\mathbf A^{\mathsf T}\mathbf A)^{-1/2}$ — ломается на вырожденных и вывернутых $\mathbf A$. Метод Müller, Bender, Chentanez, Macklin (2016) итерационный и устойчивый. На каждом шаге столбцы $\mathbf r_c$ текущего $\mathbf R$ поворачиваются к столбцам $\mathbf a_c$ матрицы $\mathbf A$ — как рамка, к которой приложен «момент»:

$$
\boldsymbol\omega = \frac{\sum_{c=1}^{3} \mathbf r_c \times \mathbf a_c}{\left|\sum_{c=1}^{3} \mathbf r_c\cdot \mathbf a_c\right| + \varepsilon}, \qquad
q \leftarrow \exp(\boldsymbol\omega)\, q .
$$

[src/math/Quaternion.h:132](../src/math/Quaternion.h#L132)
```cpp
inline Quaternion extractRotation(const Matrix3x3& A, Quaternion q, int iterations = 20) {
    for (int it = 0; it < iterations; ++it) {
        const Matrix3x3 R = q.toMatrix3x3();
        Vector3 torque(0.0f);
        float align = 0;
        for (int c = 0; c < 3; ++c) {
            torque += cross(R.col(c), A.col(c));
            align += dot(R.col(c), A.col(c));
        }
        const Vector3 omega = torque / (std::fabs(align) + 1e-9f);
        const float w = length(omega);
        if (w < 1e-9f) break;
        q = (Quaternion::fromAxisAngle(omega / w, w) * q).normalized();
    }
    return q;
}
```

> **Совет:** передавайте в `q` поворот прошлого кадра (тёплый старт). Тогда хватает нескольких итераций — `solveShapeMatching` делает 10.

---

## 1.3 Матрица 3×3 и собственные векторы (метод Якоби)

`Matrix3x3` хранится построчно (`m[row][col]`). Используется для тензоров инерции, эффективных масс контактов и поворотов. Обращение — через присоединённую матрицу; при $|\det| \le$ `singularDet` возвращается нулевая матрица (например, для двух статичных тел).

`symmetricEigen(A, eig, V)` раскладывает симметричную матрицу $\mathbf A = \mathbf V\,\operatorname{diag}(\lambda)\,\mathbf V^{\mathsf T}$ **циклическим методом Якоби** в `double`. Каждый поворот $\mathbf J_{pq}$ в плоскости $(p, q)$ обнуляет внедиагональный элемент $a_{pq}$:

$$
\theta = \frac{a_{qq} - a_{pp}}{2a_{pq}},\quad
t = \frac{\operatorname{sign}\theta}{|\theta| + \sqrt{\theta^2+1}},\quad
c = \frac{1}{\sqrt{t^2+1}},\ s = t\,c,
\qquad \mathbf A \leftarrow \mathbf J^{\mathsf T}\mathbf A\,\mathbf J,\ \ \mathbf V \leftarrow \mathbf V\mathbf J .
$$

[src/math/Matrix3x3.h:119](../src/math/Matrix3x3.h#L119)
```cpp
for (auto& pq : pairs) {
    int p = pq[0], q = pq[1];
    if (std::fabs(a[p][q]) < 1e-30) continue;
    double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
    double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
    double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
    // a <- J^T a J,  V <- V J  with J = rotation in the (p, q) plane
```

Где применяется: `ConvexHullShape` поворачивает вершины в **главные оси инерции** ([Shapes.h:107](../src/rigid/Shapes.h#L107)). Тогда тензор инерции в локальной системе тела диагонален, и хватает трёх чисел `invInertiaLocal`.

---

## 1.4 Матрица 4×4

`Matrix4x4` — аффинные преобразования и проекции. Точки — столбцы: $\mathbf p' = \mathbf M(\mathbf p, 1)$.

| Функция | Что делает |
|---|---|
| `trs(t, q, s)` | $\mathbf T\,\mathbf R\,\mathbf S$: масштаб, потом поворот, потом сдвиг ([Matrix4x4.h:33](../src/math/Matrix4x4.h#L33)) |
| `perspective(fovY, aspect, n, f)` | проекция OpenGL, камера смотрит вдоль $-z$, глубина в $[-1, 1]$ ([Matrix4x4.h:43](../src/math/Matrix4x4.h#L43)) |
| `lookAt(eye, target, up)` | матрица вида ([Matrix4x4.h:54](../src/math/Matrix4x4.h#L54)) |
| `inverse()` | общий обратный через 2×2-миноры (разложение Лапласа) ([Matrix4x4.h:96](../src/math/Matrix4x4.h#L96)) |
| `transformPoint` / `transformDirection` | с делением на $w$ / без сдвига |

---

## 1.5 Матрица N×N: LU, Холецкий, Якоби

`MatrixNxN` (в `double`) нужна для точных решений и тестов.

**LU с частичным выбором ведущего элемента.** $\mathbf P\mathbf A = \mathbf L\mathbf U$: в каждом столбце строка с наибольшим по модулю элементом становится ведущей. Множители хранятся на месте $\mathbf L$:

[src/math/MatrixNxN.cpp:71](../src/math/MatrixNxN.cpp#L71)
```cpp
for (int c = 0; c < n; ++c) {
    int piv = c; // largest entry in the column below the diagonal
    for (int i = c + 1; i < n; ++i)
        if (std::fabs(lu(i, c)) > std::fabs(lu(piv, c))) piv = i;
    if (std::fabs(lu(piv, c)) <= tiny) return false;
    if (piv != c) {
        for (int j = 0; j < n; ++j) std::swap(lu(c, j), lu(piv, j));
        std::swap(perm[c], perm[piv]);
        sign = -sign;
    }
    for (int i = c + 1; i < n; ++i) {
        lu(i, c) /= lu(c, c); // multiplier, stored in L
        for (int j = c + 1; j < n; ++j) lu(i, j) -= lu(i, c) * lu(c, j);
    }
}
```

Порог вырожденности относительный: `tiny = 1e-14 · max|a_ij|`. Затем прямой ход $\mathbf L\mathbf y = \mathbf P\mathbf b$ и обратный $\mathbf U\mathbf x = \mathbf y$ (`solveLU`, [строка 89](../src/math/MatrixNxN.cpp#L89)). На этой же базе — `determinant()` и `inverse()`.

**Холецкий** для симметричной положительно определённой матрицы ($\mathbf A = \mathbf L\mathbf L^{\mathsf T}$, [строка 136](../src/math/MatrixNxN.cpp#L136)):

$$
L_{jj} = \sqrt{a_{jj} - \sum_{k<j} L_{jk}^2}, \qquad
L_{ij} = \frac{a_{ij} - \sum_{k<j} L_{ik}L_{jk}}{L_{jj}}\ \ (i > j).
$$

Если подкоренное выражение $\le 0$, матрица не положительно определена, функция возвращает `false`.

**`solveSmall(n, M, r, x)`** ([строка 209](../src/math/MatrixNxN.cpp#L209)) — то же исключение Гаусса для $n \le 4$ во `float` и без выделения памяти. Это рабочая лошадь блочного контактного решателя (гл. 2): он решает систему до 4×4 для каждого многоточечного контакта на каждой итерации.

---

## 1.6 AABB

`AABB` — пара `lo`, `hi`. Пустой бокс имеет `lo = +∞`, `hi = −∞`, поэтому первое `expand` делает его точкой. `surfaceArea()` нужна эвристике SAH, `rayHit` — slab-тест луча:

[src/math/AABB.h:35](../src/math/AABB.h#L35)
```cpp
float rayHit(const Vector3& o, const Vector3& invDir, float tmax) const {
    float t0 = 0, t1 = tmax;
    for (int a = 0; a < 3; ++a) {
        float tn = (lo[a] - o[a]) * invDir[a];
        float tf = (hi[a] - o[a]) * invDir[a];
        if (tn > tf) std::swap(tn, tf);
        t0 = tn > t0 ? tn : t0;
        t1 = tf < t1 ? tf : t1;
        if (t0 > t1) return kInf;
    }
    return t0;
}
```

Для каждой оси луч входит в «слой» между двумя плоскостями при $t_n$ и выходит при $t_f$. Пересечение слоёв всех трёх осей — интервал $[t_0, t_1]$; если он пуст, промах.

---

## 1.7 BVH: статическая иерархия с бинированной SAH

[src/spatial/BVH.h](../src/spatial/BVH.h) — дерево, которое строится **один раз сверху вниз**. Используется для треугольных мешей (луч, ближайшая точка, знаковое расстояние), для широкой фазы CCD и в `BvhBroadPhase`.

### Эвристика площади поверхности (SAH)

Вероятность того, что случайный луч попадёт в дочерний бокс, пропорциональна его площади. Поэтому стоимость разбиения узла на левую ($L$) и правую ($R$) части

$$
C = N_L\,S(L) + N_R\,S(R),
$$

где $N$ — число примитивов, $S$ — площадь поверхности бокса. Перебирать все возможные разрезы дорого, поэтому центроиды раскладываются по **12 корзинам** вдоль каждой оси, и оцениваются 11 разрезов между корзинами за один проход вперёд и один назад:

[src/spatial/BVH.cpp:57](../src/spatial/BVH.cpp#L57)
```cpp
for (int b = 0; b < kBins - 1; ++b) {
    acc.expand(binBox[b]);
    n += binCount[b];
    leftArea[b] = acc.valid() ? acc.surfaceArea() : 0.0f;
    leftCount[b] = n;
}
acc = AABB();
n = 0;
for (int b = kBins - 1; b > 0; --b) {
    acc.expand(binBox[b]);
    n += binCount[b];
    float cost = leftCount[b - 1] * leftArea[b - 1] + n * (acc.valid() ? acc.surfaceArea() : 0.0f);
    if (leftCount[b - 1] > 0 && n > 0 && cost < bestCost) {
```

Узел становится листом, если примитивов ≤ `maxLeafSize` или если разрез не дешевле листа ($C \ge N\,S$) при $N \le 16$.

### Знаковое расстояние до меша: псевдонормали с угловыми весами

`MeshBVH::closestPoint` находит ближайшую точку на меше обходом дерева (ближний потомок первым, отсечение по `distance2` до бокса). Сложнее определить **знак**: внутри точка или снаружи. Нормаль ближайшего треугольника врёт, если ближайшая точка лежит на ребре или вершине.

Решение — псевдонормали (Bærentzen & Aanæs 2005). У каждого элемента своя нормаль:

- грань — её нормаль;
- ребро — сумма нормалей двух смежных граней;
- вершина — сумма нормалей смежных граней, **взвешенная углом** грани при этой вершине.

С такими нормалями знак $\operatorname{sign}\big((\mathbf p - \mathbf p_{closest})\cdot\mathbf n_{pseudo}\big)$ точен для замкнутого меша. Функция `closestPtTri` (Ericson, *Real-Time Collision Detection*, §5.1.5) возвращает, какой элемент ближайший: грань, одна из вершин или одно из рёбер.

[src/spatial/BVH.cpp:262](../src/spatial/BVH.cpp#L262)
```cpp
Vector3 N;
const auto& tri = tris_[bestTri];
switch (bestRegion) {
case 0: N = faceN_[bestTri]; break;
case 1: case 2: case 3: N = vertexN_[tri[bestRegion - 1]]; break;
default: N = edgeN_[size_t(bestTri) * 3 + (bestRegion - 4)]; break;
}
```

---

## 1.8 Динамическое AABB-дерево

[src/spatial/AABBTree.h](../src/spatial/AABBTree.h) — иерархия, которая **обновляется по ходу движения** (как `b2DynamicTree` в Box2D и `btDbvt` в Bullet; E. Catto, *Dynamic Bounding Volume Hierarchies*, GDC 2019).

```mermaid
flowchart LR
    A["update(proxy, box)"] --> B{"box внутри<br/>«толстого» бокса?"}
    B -- да --> C["ничего не делать"]
    B -- нет --> D["removeLeaf"] --> E["толстый бокс = box + запас"] --> F["insertLeaf (SAH)"] --> G["refitUpwards:<br/>боксы, высоты, AVL-повороты"]
```

### Толстые боксы

Лист хранит бокс объекта, расширенный на `fatten × размер` с каждой стороны ([AABBTree.cpp:18](../src/spatial/AABBTree.cpp#L18)). Пока объект дрожит на месте (стопка, куча), его настоящий бокс остаётся внутри толстого, и дерево не меняется вовсе. В тесте «jiggling in place» из 672 объектов перевставлено только 7.

### Вставка: спуск по SAH

Новый лист `L` нужно сделать братом какого-то узла `S` (под ними появится новый родитель). Спускаемся от корня и в каждом узле сравниваем:

- **создать родителя здесь**: стоимость $2\,S(N \cup L)$;
- **спуститься в ребёнка $C$**: каждый предок и так вырастет на $\Delta = 2\,(S(N\cup L) - S(N))$ («унаследованная» стоимость), плюс рост самого ребёнка.

[src/spatial/AABBTree.cpp:101](../src/spatial/AABBTree.cpp#L101)
```cpp
while (!nodes_[index].isLeaf()) {
    const Node& n = nodes_[index];
    const float area = n.box.surfaceArea();
    const float combined = unite(n.box, box).surfaceArea();
    // Making a new parent for this node and the leaf here costs:
    const float costHere = 2.0f * combined;
    // Going further down, every node on the way grows by this much:
    const float inherited = 2.0f * (combined - area);
    auto costDown = [&](int child) {
        const Node& c = nodes_[child];
        float grown = unite(box, c.box).surfaceArea();
        return (c.isLeaf() ? grown : grown - c.box.surfaceArea()) + inherited;
    };
    const float costLeft = costDown(n.left), costRight = costDown(n.right);
    if (costHere < costLeft && costHere < costRight) break;
    index = costLeft < costRight ? n.left : n.right;
}
```

Удаление проще: место родителя удаляемого листа занимает его брат ([AABBTree.cpp:140](../src/spatial/AABBTree.cpp#L140)).

### Балансировка поворотами (как в AVL-дереве)

После каждой вставки или удаления путь до корня «переподгоняется» (`refitUpwards`), и в каждом узле вызывается `balance`. Если один ребёнок выше другого больше чем на 1 уровень, он поворачивается вверх:

```
          A                    C
        /   \                /   \
       B     C      ->      A     F     (если G выше F — наоборот)
            / \            / \
           F   G          B   G
```

[src/spatial/AABBTree.cpp:179](../src/spatial/AABBTree.cpp#L179) — код поворота. Высокий ребёнок `U` занимает место `A`; `A` оставляет себе второго ребёнка и забирает у `U` его **низкого** внука. Так высота дерева остаётся $O(\log n)$. В тесте 672 объекта дают высоту 11 при $\log_2 n = 9.4$.

### Запросы

`query(box, fn)` и `raycast(o, d, maxT, fn)` — обход стеком с отсечением по боксу. `findPairs` собирает все пары с пересекающимися толстыми боксами.

---

## Проверка

| Тест (`RF_TEST=...`) | Что проверяет |
|---|---|
| `math` | $\mathbf A\mathbf A^{-1} = \mathbf I$, восстановление $\mathbf V\Lambda\mathbf V^{\mathsf T}$ с ошибкой < 1e-5; $\exp(\log q) = q$; slerp на полпути даёт половину угла; LU и Холецкий на матрице типа Гильберта 6×6 с ошибкой < 1e-12; невязка собственных векторов N×N < 1e-10; `solveSmall` 4×4 |
| `mass properties` | объём и инерция бокса как многогранника; Якоби диагонализует повёрнутый тензор |
| `bvh` | 2000 случайных точек: ближайшая точка совпадает с перебором, знак расстояния совпадает с уравнением эллипсоида (0 ошибок); луч по сфере, $t = 4$ |
| `dynamic AABB tree` | инварианты дерева (`validate()`), запросы = перебор, высота ≤ $2\log_2 n + 2$, дрожащие объекты не трогают дерево |
| `broad phase` | BVH, SAP и AABB-дерево дают **ровно те же пары**, что перебор, 200 кадров подряд |

## Литература

- M. Müller, J. Bender, N. Chentanez, M. Macklin. *A Robust Method to Extract the Rotational Part of Deformations.* MIG 2016.
- C. G. J. Jacobi (1846); G. H. Golub, C. F. Van Loan. *Matrix Computations*, 4th ed., 2013 (§8.5 — метод Якоби).
- I. Wald. *On fast Construction of SAH-based Bounding Volume Hierarchies.* IEEE RT 2007 (бинированная SAH).
- J. A. Bærentzen, H. Aanæs. *Signed Distance Computation Using the Angle Weighted Pseudonormal.* IEEE TVCG 11(3), 2005.
- C. Ericson. *Real-Time Collision Detection.* Morgan Kaufmann, 2005.
- E. Catto. *Dynamic Bounding Volume Hierarchies.* GDC 2019.
