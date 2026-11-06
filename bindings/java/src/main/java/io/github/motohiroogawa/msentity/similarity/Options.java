package io.github.motohiroogawa.msentity.similarity;

import static io.github.motohiroogawa.msentity.similarity.Native.call;
import static io.github.motohiroogawa.msentity.similarity.Native.check;

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Set;

/**
 * Score settings. Immutable; change values with the {@code with…} methods:
 *
 * <pre>{@code
 * Options options =
 *     Options.defaults()
 *         .withMethods(Method.COSINE, Method.MODIFIED_COSINE)
 *         .withMzTolerancePpm(10);
 * }</pre>
 *
 * @param methods the scores to compute (at least one); result tables get one column per method,
 *     named {@code cosine}, {@code reverse_cosine}, {@code modified_cosine}, {@code bonanza} (in
 *     this order); the precursor methods give null (NaN from pairScores / allPairs) when a
 *     spectrum has no precursor m/z. {@link Similarity#pairScores} and {@link Similarity#allPairs}
 *     need exactly one.
 * @param mzToleranceDa absolute m/z tolerance in Da; 0 disables it
 * @param mzTolerancePpm symmetric m/z tolerance in ppm; 0 disables it
 * @param intensityExponent intensities are raised to this power (&gt; 0) before scoring
 * @param threads worker threads; 0 = all hardware threads
 * @param precursorColumn precursor m/z column for the precursor methods; null or empty detects
 *     PrecursorMZ, precursor_mz, PEPMASS, or ParentMZ
 */
public record Options(Set<Method> methods, double mzToleranceDa, double mzTolerancePpm,
                      double intensityExponent, int threads, String precursorColumn) {
  /** A scoring method. */
  public enum Method {
    /** Tolerance-based cosine similarity. */
    COSINE(0x1),
    /** Cosine where the query norm counts only peaks matched in the reference. */
    REVERSE_COSINE(0x2),
    /** Cosine that also pairs peaks shifted by the precursor m/z difference; needs precursors. */
    MODIFIED_COSINE(0x4),
    /** Bonanza score over direct and precursor-shifted pairs; needs precursors. */
    BONANZA(0x8);

    final int flag;

    Method(int flag) {
      this.flag = flag;
    }
  }

  public Options {
    if (methods == null || methods.isEmpty()) {
      throw new SimilarityException.InvalidArgument("at least one method is required");
    }
    methods = Collections.unmodifiableSet(EnumSet.copyOf(methods));
  }

  /** cosine, 0.01 Da or 10 ppm tolerance, exponent 1, all threads, detected precursor column. */
  public static Options defaults() {
    return new Options(EnumSet.of(Method.COSINE), 0.01, 10.0, 1.0, 0, null);
  }

  public Options withMethods(Method first, Method... rest) {
    return withMethods(EnumSet.of(first, rest));
  }

  public Options withMethods(Set<Method> value) {
    return new Options(value, mzToleranceDa, mzTolerancePpm, intensityExponent, threads,
                       precursorColumn);
  }

  /** Methods from comma-separated names, e.g. {@code "cosine,bonanza"}. */
  public Options withMethods(String text) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment options = toNative(arena);
      MemorySegment t = Native.string(arena, text);
      check(call(() -> (int) Native.parseMethods.invoke(t, options)));
      int flags = options.get(ValueLayout.JAVA_INT, 0);
      EnumSet<Method> parsed = EnumSet.noneOf(Method.class);
      Arrays.stream(Method.values()).filter(m -> (flags & m.flag) != 0).forEach(parsed::add);
      return withMethods(parsed);
    }
  }

  public Options withMzToleranceDa(double value) {
    return new Options(methods, value, mzTolerancePpm, intensityExponent, threads, precursorColumn);
  }

  public Options withMzTolerancePpm(double value) {
    return new Options(methods, mzToleranceDa, value, intensityExponent, threads, precursorColumn);
  }

  /**
   * Both tolerances from text such as {@code "0.01Da,10ppm"} or {@code "10ppm"}: {@code <value>Da}
   * and/or {@code <value>ppm} in either order, units case-insensitive; a unit left out becomes 0.
   */
  public Options withMzTolerance(String text) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment options = toNative(arena);
      MemorySegment t = Native.string(arena, text);
      check(call(() -> (int) Native.parseMzTolerance.invoke(t, options)));
      return new Options(methods, options.get(ValueLayout.JAVA_DOUBLE, 8),
                         options.get(ValueLayout.JAVA_DOUBLE, 16), intensityExponent, threads,
                         precursorColumn);
    }
  }

  /** Both tolerances as canonical text, e.g. {@code "0.01Da,10ppm"} (a zero term is left out). */
  public String formatMzTolerance() {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment options = toNative(arena);
      MemorySegment slot = arena.allocate(ValueLayout.ADDRESS);
      check(call(() -> (int) Native.formatMzTolerance.invoke(options, slot)));
      return Native.takeString(slot);
    }
  }

  public Options withIntensityExponent(double value) {
    return new Options(methods, mzToleranceDa, mzTolerancePpm, value, threads, precursorColumn);
  }

  public Options withThreads(int value) {
    return new Options(methods, mzToleranceDa, mzTolerancePpm, intensityExponent, value,
                       precursorColumn);
  }

  public Options withPrecursorColumn(String value) {
    return new Options(methods, mzToleranceDa, mzTolerancePpm, intensityExponent, threads, value);
  }

  /** struct mssim_options in {@code arena}. */
  MemorySegment toNative(Arena arena) {
    int flags = 0;
    for (Method method : methods) flags |= method.flag;
    MemorySegment out = arena.allocate(Native.OPTIONS);
    out.set(ValueLayout.JAVA_INT, 0, flags);
    out.set(ValueLayout.JAVA_INT, 4, threads);
    out.set(ValueLayout.JAVA_DOUBLE, 8, mzToleranceDa);
    out.set(ValueLayout.JAVA_DOUBLE, 16, mzTolerancePpm);
    out.set(ValueLayout.JAVA_DOUBLE, 24, intensityExponent);
    out.set(ValueLayout.ADDRESS, 32, Native.string(arena, precursorColumn));
    return out;
  }
}
