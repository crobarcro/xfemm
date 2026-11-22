# XFEMM licensing

## Project licence

Except for files and components that carry a separate third-party licence, the XFEMM code in this repository is offered under the **GNU General Public License, version 3, or (at your option) any later version** (`GPL-3.0-or-later`). The complete GPL version 3 text is in [`LICENSE`](LICENSE).

This GPL grant covers the FEMM-derived code for which David C. Meeker supplied the additional licence grant reproduced below, together with XFEMM contributions for which the relevant copyright holders have authorised relicensing. The XFEMM maintainer has obtained the required relicensing permissions from the identified XFEMM copyright holders/contributors.

Some source files still contain historical Aladdin Free Public License (AFPL) notices. Those notices document the earlier licensing history and do not negate the additional GPL permissions granted by the relevant copyright holders. In particular, David C. Meeker's grant expressly leaves the AFPL in place as an alternative for his contributions; it does not revoke or modify copies previously distributed under the AFPL.

## David C. Meeker's additional licence grant

The following grant is reproduced verbatim:

> **GRANT OF ADDITIONAL LICENSE**
>
> I, David C. Meeker of Natick, Massachusetts, am the author and copyright holder of the Finite Element Method Magnetics (FEMM) source code from which portions of the xfemm project are derived.
>
> I hereby grant that my contributions to FEMM, as incorporated into the xfemm project, may additionally be used, copied, modified, and distributed under the terms of the GNU General Public License, version 3, or (at the recipient's option) any later version. This grant is non-exclusive, irrevocable, and royalty-free, and extends to all recipients of xfemm and of works derived from it.
>
> Relationship to the Aladdin Free Public License: this grant is made in addition to, and not in substitution for, the AFPL under which this code has previously been distributed. Nothing here revokes, terminates, or modifies that license, and copies already distributed under it are unaffected. A recipient who accepts the terms of this grant may rely on the GNU General Public License alone; such a recipient has no obligation under the Aladdin Free Public License with respect to my contributions.

## Third-party components

The GPL relicensing does **not** relicense third-party material for which the XFEMM contributors do not control copyright. A file or component with its own licence notice remains subject to that notice. Important examples include:

- **Triangle**, including the bundled Triangle sources used by the Triangle mesher backend: see [`cfemm/LICENSE-triangle.txt`](cfemm/LICENSE-triangle.txt). Triangle is not relicensed under the GPL by this change.
- **Lua** sources bundled with `cfemm`: see [`cfemm/LICENSE-Lua.txt`](cfemm/LICENSE-Lua.txt).
- Other files that contain an explicit separate third-party licence or copyright notice, including third-party helper code under `mfemm`, retain those terms.

The historical [`cfemm/LICENSE-FEMM.txt`](cfemm/LICENSE-FEMM.txt) and `cfemm/fmesher/LICENCE.txt` files are retained to document the pre-existing FEMM/AFPL licensing terms. For the additional GPL licence now offered for eligible XFEMM and FEMM-derived code, use the repository-level [`LICENSE`](LICENSE) together with this file.
