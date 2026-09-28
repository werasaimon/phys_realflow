# PhysRealFlow

Физический SDK на **C++17 без внешних зависимостей**; редактор на Qt 6 — в [отдельном репозитории](https://github.com/werasaimon/phys_realflow_editor).
Всё в единицах СИ, каждый метод — со ссылкой на статью и с численным тестом.

| Подсистема | Методы |
|---|---|
| Твёрдые тела | GJK/EPA, SAT, блочный LCP, CCD, сочленения, выпуклая декомпозиция |
| Частицы (как NVIDIA FleX) | жидкость PBF, мягкие тела, ткань XPBD с разрывом по нитям и швам |
| Газ | Навье–Стокс на MAC-сетке, аэротруба, дым, двусторонняя связь с телами |
| Огонь | горение с кислородным пределом, излучение пламени, теплопроводность, пиролиз хлопка по Аррениусу |
| Плазма | резистивная МГД: сила Лоренца, закон индукции (constrained transport), магнитосфера |
| Рендер | OpenGL 3.0 (работает и без видеокарты): пламя как чёрное тело, экранная вода, светящаяся плазма |

**Документация** — физика, численные методы, код и проверка по главам: [docs/README.md](docs/README.md).
**Замысел** — для кого это, мерка качества (граница Крамера–Рао) и карта пути: [docs/00-vision.md](docs/00-vision.md).

## Быстрый старт

```sh
# SDK и тесты (любой компилятор C++17, CMake ≥ 3.21, никаких зависимостей)
cmake -S . -B build-core -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-core -j
build-core/rf_tests
```

Редактор с окном, панелями и графиками — отдельный репозиторий [phys_realflow_editor](https://github.com/werasaimon/phys_realflow_editor),
в котором этот SDK подключён git-сабмодулем (`git clone --recurse-submodules`).

## Структура

```
src/      SDK (rfcore) — только стандартная библиотека C++:
            math/ core/ spatial/   векторы, кватернионы, меши, пул потоков, BVH
            rigid/                 RigidWorld + ContactSolver, Islands, ShockPropagation, Grab, Joints, TimeOfImpact
            particles/             ParticleSystem + DensitySolver, ParticleContacts; SoftBody, Cloth
            gas/                   GasSolver + Advection, PressureSolver, MovingSolids, Heat; Combustion
            plasma/                MagneticField
            scene/                 Simulation (фасад) + Obstacle, Coupling, Snapshot; интерфейс Scene
            core/Probe             отладчик: любой канал по имени (значения, счётчики, таймеры), отладочная отрисовка
samples/  готовые сцены (rfsamples) — решения на движке, как samples/ у Box2D:
            Liquid, WindTunnel, Smoke, Rigid (+ рельеф), SoftBody, Fire, Hydro; plasma/ Tokamak, Magnetosphere
tests/    rf_tests — численные эксперименты с аналитическим ответом
docs/     документация с формулами, кодом и графиками
```

## Лицензия

PhysRealFlow — общее достояние. Лицензия MIT-0 ([LICENSE](LICENSE)): берите и делайте что хотите —
изучайте, меняйте, встраивайте, продавайте, — даже упоминать нас не обязательно. Единственное
исключение — чужой код, перенесённый из других проектов: он остаётся под своими лицензиями, и они
перечислены в [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Поддержать проект

Инструмент бесплатный и таким останется. Если он вам помог — в учёбе, в исследовании или в работе, —
можно по желанию поддержать его разработку пожертвованием. Это добровольно и ни на что не влияет:
все возможности доступны всем.
