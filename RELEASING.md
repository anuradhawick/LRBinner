# Releasing LRBinner

This document describes how to publish LRBinner to [TestPyPI](https://test.pypi.org) and [PyPI](https://pypi.org) using the GitHub Actions workflow in `.github/workflows/publish.yml`.

---

## 1. One-time setup: link PyPI environments via OIDC trusted publishing

The workflow uses [PyPA's OIDC trusted publishing](https://docs.pypi.org/trusted-publishers/) — no API tokens required. You must configure two environments in both GitHub and PyPI/TestPyPI.

### 1a. Create GitHub environments

1. Go to **Settings → Environments** in the GitHub repository.
2. Click **New environment** and create two environments:

   | Environment name | Purpose |
   |---|---|
   | `pypi` | Production releases to [pypi.org](https://pypi.org) |
   | `testpypi` | Pre-release testing on [test.pypi.org](https://test.pypi.org) |

3. For the `pypi` environment, add a **required reviewer** (manual approval gate):
   - Under **Environment protection rules**, enable **Required reviewers**.
   - Add yourself (or a team) as a required reviewer.
   - This means every push to PyPI requires an explicit approval click before the job runs.

   The `testpypi` environment does **not** need required reviewers — it is triggered manually via `workflow_dispatch` and is used for testing only.

### 1b. Configure a trusted publisher on TestPyPI

1. Log in to [https://test.pypi.org](https://test.pypi.org).
2. Go to **Account settings → Publishing → Add a new publisher**.
3. Fill in:

   | Field | Value |
   |---|---|
   | PyPI project name | `LRBinner` |
   | Owner | `anuradhawick` |
   | Repository name | `LRBinner` |
   | Workflow name | `publish.yml` |
   | Environment name | `testpypi` |

4. Click **Add**.

### 1c. Configure a trusted publisher on PyPI

1. Log in to [https://pypi.org](https://pypi.org).
2. Go to **Account settings → Publishing → Add a new publisher**.
3. Fill in:

   | Field | Value |
   |---|---|
   | PyPI project name | `LRBinner` |
   | Owner | `anuradhawick` |
   | Repository name | `LRBinner` |
   | Workflow name | `publish.yml` |
   | Environment name | `pypi` |

4. Click **Add**.

> **First release only:** If the project doesn't yet exist on PyPI/TestPyPI you must use *pending publisher* mode — select "PyPI project doesn't exist yet" when adding the publisher. PyPI will create the project automatically on the first successful publish.

---

## 2. Releasing to TestPyPI (manual dispatch)

TestPyPI releases are triggered manually at any time — useful for verifying the package metadata, install flow, and wheel before a real release.

1. Go to **Actions → Publish to PyPI** in the GitHub repository.
2. Click **Run workflow** (top-right of the workflow list).
3. Select the branch you want to publish from (e.g. `main`) and click **Run workflow**.
4. The workflow will:
   - Build the wheel and sdist with `flit build`.
   - Publish to TestPyPI under the `testpypi` environment (no approval required).
5. Verify the release at `https://test.pypi.org/project/LRBinner/`.
6. Test the install from TestPyPI:

   ```bash
   pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ LRBinner
   ```

---

## 3. Releasing to PyPI (tag-based, manual approval)

Production releases are triggered by pushing a version tag. The `pypi` environment's required-reviewer rule means the publish job is gated behind a manual approval click.

### Step 1 — bump the version

Update `__version__` in `lrbinner/__init__.py`:

```python
__version__ = "2.2.0"   # new version
```

Commit the change:

```bash
git add lrbinner/__init__.py
git commit -m "chore: bump version to 2.2.0"
git push
```

### Step 2 — create and push a version tag

```bash
git tag v2.2.0
git push origin v2.2.0
```

The tag must match the pattern `v*` (e.g. `v2.2.0`, `v3.0.0-rc1`).

### Step 3 — approve the PyPI deployment

1. After the tag is pushed, GitHub starts the **Publish to PyPI** workflow automatically.
2. The `build` job runs immediately and uploads the wheel/sdist as a workflow artifact.
3. The `publish-pypi` job then **pauses** and waits for a required reviewer to approve.
4. Go to **Actions → Publish to PyPI → (the running workflow)**.
5. Click **Review deployments**, select `pypi`, and click **Approve and deploy**.
6. The job uploads the artifacts to PyPI using OIDC — no token exchange needed.

### Step 4 — verify

Check the release at `https://pypi.org/project/LRBinner/` and verify the install:

```bash
pip install LRBinner==2.2.0
lrbinner --version
```

---

## Workflow summary

```
push tag v*
    └── build job (flit build → artifact)
            ├── publish-pypi  [requires manual approval via 'pypi' environment]
            └── (publish-testpypi skipped — only runs on workflow_dispatch)

workflow_dispatch (manual)
    └── build job (flit build → artifact)
            ├── publish-testpypi  [runs immediately, no approval]
            └── (publish-pypi skipped — only runs on tag push)
```
