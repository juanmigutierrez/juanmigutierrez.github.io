# Deploying the interactive blog post to juanmigutierrez.github.io

The page is **fully self-contained** (real 3D data is embedded; three.js loads from a CDN).
It matches your site's existing design system (Syne / DM Mono / DM Sans, dark theme).

## 1. Copy two things into your `juanmigutierrez.github.io` repo

```
superquadric-fitting.html   ->  repo root  (next to index.html)
assets/img/*.png            ->  assets/img/   (8 result images)
```

The HTML references images as `assets/img/...`, so keeping that folder structure "just works"
alongside your existing `assets/profile.png`.

## 2. Add a card on your homepage

Paste this into `index.html`, e.g. as the first card inside the `04 — Selected Projects` list:

```html
<a class="project-card" href="superquadric-fitting.html">
  <div class="project-header">
    <span class="project-name">Robust 3D Primitive Fitting — interactive</span>
    <span class="project-lang">Three.js</span>
  </div>
  <p class="project-desc">Drag real fitted point clouds in 3D. Extending GAIR-RANSAC (STAG 2025) to recover superquadrics from noisy, sharp, overlapping data.</p>
</a>
```

## 3. Commit & push

```bash
git add superquadric-fitting.html assets/img
git commit -m "Add interactive superquadric-fitting write-up"
git push
```

Live at: https://juanmigutierrez.github.io/superquadric-fitting.html

## Notes
- The 3D viewer needs WebGL + internet (three.js r160 via jsDelivr). A graceful fallback
  message appears if unavailable; the static result images below still tell the story.
- To embed inside another page instead, wrap it in an `<iframe src="superquadric-fitting.html">`.
