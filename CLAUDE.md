# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

PhysRealFlow — физический SDK на C++17: твёрдые тела, частицы (жидкость, мягкие тела, ткань),
газ на MAC-сетке, огонь, резистивная МГД, геодезические в метрике Керра. Всё в СИ. Qt-редактор
живёт в отдельном репозитории `phys_realflow_editor` и подключает этот SDK сабмодулем, поэтому
здесь нет ни Qt, ни рендера.

## Сборка и тесты

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
build/rf_tests                              # код выхода = число проваленных CHECK
RF_TEST="box stack" build/rf_tests          # только тесты, в имени которых есть подстрока
```

Имена тестов и порядок запуска — в `tests/main.cpp` (`run("имя", fn)`); новый тест добавляется
туда и в список `rf_tests` в `CMakeLists.txt`. Отдельного фреймворка нет: `CHECK(cond, fmt, ...)`
из `tests/TestRunner.h`.

Как в CI (GitHub Actions и `.gitlab-ci.yml`, одинаковые проверки):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRF_WERROR=ON -DRF_STRICT_FP=ON
# санитайзеры: -DCMAKE_BUILD_TYPE=RelWithDebInfo -DRF_SANITIZE=ON -DRF_STRICT_FP=ON,
# гоняются только быстрые тесты (список в .github/workflows/ci.yml), ASAN_OPTIONS=detect_leaks=0
```

Опции CMake: `RF_WERROR`, `RF_STRICT_FP`, `RF_SANITIZE`, `RF_LTO`, `RF_BUILD_TESTS`, `RF_BUILD_VERIFY`.

Переменные окружения:

| Переменная | Что делает |
|---|---|
| `RF_TEST` | фильтр по подстроке имени; также включает тесты «только по запросу» (tokamak, `determinism: the fire`, `benchmark: gas pressure`) |
| `RF_THREADS` | размер пула потоков (`src/core/Parallel.cpp`) |
| `RF_PLOT_DIR` | куда тесты пишут CSV для графиков документации |
| `RF_JUNIT` | JUnit XML для CI |
| `RF_LONG_ORBITS`, `RF_TOKAMAK_RES`, `RF_WRITE_IMAGES` | параметры отдельных длинных тестов |

Проверка читаемости (обязана проходить): `python3 tools/readability.py .` — функция не длиннее
60 строк, каждый файл в `src/`, `samples/`, `verification/` начинается с абзаца-комментария.

## Верификация (`rf_verify`, ASME V&V 20)

`verification/` — реестр эталонных случаев (MMS, сходимость по Ричардсону, метрика валидации,
Монте-Карло по входам), цель `rf_verify`. Подробно — `docs/10-verification.md`.

```sh
build/rf_verify --quick | --full | --case <id> | --list | --replot verification/results/<file>.json
```

- Каждый прогон пишет `verification/results/<UTC>-<commit>.json`; `--full` и `--board` перегенерируют
  `docs/09-benchmarks.md` (только сгенерированный блок) и `docs/img/verification/`.
- Слепой анализ: у части случаев эталон запечатан (`verification/sealed/`). `--seal`, `--freeze`,
  `--unblind` и журнал `verification/unblinding-log.md` — процедура ответственного за релиз;
  не запускать их и не подгонять решатель под эталон без прямой просьбы.

## Архитектура

Зависимости только «вниз»: `math` → `core`/`spatial` → решатели (`rigid`, `particles`, `gas`,
`plasma`, `relativity`) → `scene` (фасад). Три CMake-цели:

- **`rfcore`** (`src/`) — только стандартная библиотека C++17 и `Threads::Threads`. Конфигурация
  CMake падает, если к `rfcore` подключено что-то ещё. Файлы добавляются в `CMakeLists.txt` вручную.
- **`rfsamples`** (`samples/`) — готовые сцены, как `samples/` у Box2D. Ядро не знает ни одной сцены.
  Номер сцены — значение `enum class Preset` в `samples/Samples.h`; порядок не менять, на номера
  ссылаются документация и командная строка редактора (`--preset N`).
- **`rfverify`** + `rf_verify` (`verification/`) — реестр V&V поверх `rfcore`.

Ключевые места, которые требуют чтения нескольких файлов:

- **Фасад `Simulation`** (`src/scene/Simulation.h`) держит три решателя (`rigid`, `particles`, `grid`),
  препятствие, шаг кадра со связками между решателями (`Coupling.cpp`) и `RenderSnapshot` (`Snapshot.cpp`).
  Сцена (`src/scene/Scene.h`) получает его в `configure` → `build` (при каждом `reset`) → `afterStep` → `describe`.
- **Большой класс разрезан по ответственностям**, как в Box2D: методы `RigidWorld` лежат в
  `RigidWorld.cpp`, `ContactSolver.cpp`, `Islands.cpp`, `ShockPropagation.cpp`, `Grab.cpp`, `XpbdSolver.cpp`…;
  `GasSolver` — в `Advection.cpp`, `PressureSolver.cpp`, `Multigrid.cpp`, `MovingSolids.cpp`, `Heat.cpp`…;
  `ParticleSystem` — в `DensitySolver.cpp`, `ParticleContacts.cpp`. Метод ищется по ответственности, а не по имени класса.
- **Граф сцены** (`SceneGraph`, `GraphScene`, `MetaObjects.cpp`, `EntityShapes.cpp`) — сцена как данные,
  так её собирает редактор. Три слоя: источник (геометрия) → сущность (поза + роли) → мета-объекты
  в решателях. `rebuildEntity` всегда строит заново из источника, а не из текущего состояния.
  Из ролей «из чего сделано» (`rigid`/`soft`/`liquid`/`cloth`) действует одна; `magnet`, `emitter`,
  `flammable`, `heat` комбинируются. Текстовый формат графа обязан давать `save → load → save` до бита.
- **`Probe`** (`src/core/Probe.h`) — глобальный канал отладки: значения, счётчики, таймеры и отладочная
  отрисовка по слоям. Имена каналов `"part/quantity"` в ASCII; писать раз в кадр или подшаг, не на частицу.
  Тестовая программа подменяет `operator new`, чтобы считать выделения за кадр; SDK сам аллокатор не трогает.

## Инварианты, которые ловят тесты

- **Детерминизм.** Два прогона и прогоны на 1, 2 и всех потоках дают одинаковые биты. С `RF_STRICT_FP=ON`
  хэш эталонной сцены твёрдых тел обязан совпасть с `kGolden` в `tests/DeterminismTests.cpp`.
  Константа снята на MinGW (GCC 11.2): `std::sin` в glibc и MSVC округляет последний бит иначе,
  поэтому на Linux расхождение хэша само по себе не значит, что решатель сломан. Если результаты
  твёрдых тел меняются намеренно, обновляется `kGolden`, дописывается строка в историю хэшей над
  тестом, а причина объясняется в коммите.
- Тригонометрия через `rf::sin`, `rf::cos`, `rf::atan2` (`src/math/ElementaryFunctions.h`); тест
  сверяет их со `std::` до бита.

## Правила кода и документации

- Один класс — один файл с тем же именем. Никакой шаблонной магии: обычные классы, `std::vector`, `std::function`.
- **Код не зависит от платформы.** Только стандартный C++17, без `#ifdef _WIN32` и системных API.
  Пишется на Windows (MinGW), собираться и работать обязан так же на Linux и macOS.
- **Числа в файлах** читаются и пишутся только через `rf::numberText` и `rf::parseNumber` из
  `src/core/Format.h`, а не через `printf` и `strtof`. Qt на Linux включает системную локаль, и в
  русской `printf` пишет `9,81`, а `strtof` не читает `9.81`. `rf::format` годится только для показаний людям.
  Это держит тест `scene graph: the file does not depend on the locale`.
- Каждая формула в коде — с комментарием и ссылкой на статью; каждый метод — с численным тестом
  против аналитического ответа или эталона.
- Документация в `docs/` на русском: физика → метод → код со ссылкой `файл#Lстрока` → тест.
  Графики не рисуются руками: тест пишет CSV в `RF_PLOT_DIR`, скрипт `tools/plot_*.py` (только
  стандартная библиотека Python) делает SVG. Команда для каждого графика — в docstring скрипта.
  При сдвиге кода ссылки `#L<строка>` в `docs/` надо проверить.
- Код комментируется по-английски, коммиты и документация — по-русски.
- Чужой код, перенесённый из других проектов, вносится в `THIRD_PARTY_NOTICES.md` с его лицензией.
